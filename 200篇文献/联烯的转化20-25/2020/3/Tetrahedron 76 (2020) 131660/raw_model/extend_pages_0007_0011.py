#!/usr/bin/env python3
"""Evidence-linked v3.1 extraction for PDF pages 7-11.

The script replaces only the existing records whose scheme summaries are
expanded by experimental procedures on pages 7-10, and augments the HWE
record with the page-11 footnote-39 side product.  It emits a merge bundle;
it does not edit schemas, validators, or chemistry rules.
"""
from copy import deepcopy
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_c724145241bb10a8"
DATASET = st.read(RUN / "dataset.json")
OLD = st.read(RUN / "standard_records.json")
OLD_M = {m["molecule_id"]: m for m in OLD["molecules"]}
OLD_R = {r["reaction_id"]: r for r in OLD["reactions"]}
E = {p: [f"{DOC}:page{p}:image", f"{DOC}:page{p}:text"] for p in range(1, 12)}


def ev(*pages):
    return list(dict.fromkeys(x for p in pages for x in E[p]))


def F(raw=None, norm=None, evidence=(), method="original_extraction", status=None,
      verification=None, derivation=None, alternatives=()):
    return st.fact(raw, norm, evidence, method=method, status=status,
                   verification=verification, derivation=derivation,
                   alternatives=alternatives)


def Q(number, unit, raw=None, qualifier=None):
    v = {"number": number, "unit": unit, "raw_text": raw or f"{number} {unit}"}
    if qualifier:
        v["qualifier"] = qualifier
    return raw or v["raw_text"], v


def Rng(lo, hi, unit, raw):
    return raw, {"minimum": lo, "maximum": hi, "unit": unit, "raw_text": raw}


def P(pid, mid, role, evidence, stage="reaction", **charges):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    for key, val in charges.items():
        if key == "reference" or val is None:
            continue
        raw, norm = val
        p["charges"][key] = F(raw, norm, evidence)
    if charges.get("equivalents") is not None or charges.get("mol_percent") is not None:
        ref = charges.get("reference")
        p["equivalence_reference"] = F(ref, ref, evidence)
    return p


def S(sid, evidence, pids, operation, temperature=None, time=None, stage="reaction",
      atmosphere=None):
    s = st.new_step(sid, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = F(operation, operation, evidence)
    if temperature is not None:
        raw, norm = temperature if isinstance(temperature, tuple) else (temperature, temperature)
        s["temperature"] = F(raw, norm, evidence)
    if time is not None:
        raw, norm = time if isinstance(time, tuple) else (time, time)
        s["time"] = F(raw, norm, evidence)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, evidence)
    return s


def O(oid, kind, targets, evidence, value=None, basis=None, yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    if value is not None:
        raw, norm = value if isinstance(value, tuple) else (value, value)
        o["value"] = F(raw, norm, evidence)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, evidence)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, evidence)
    if scope is not None:
        o["scope"] = F(scope, scope, evidence)
    return o


def rxn(suffix, pages, participants, steps, outcomes, transformation, unresolved=(),
        notes=(), location=None, scale=None, solvent=None, classification=None):
    rid = f"{DOC}:{suffix}"
    r = deepcopy(OLD_R[rid])
    evidence = ev(*pages)
    r["evidence_ids"] = evidence
    r["location"] = F(location or f"PDF pages {', '.join(map(str, pages))}",
                        {"pages": list(pages), "region": location or suffix}, evidence)
    r["transformation_type"] = F(transformation, transformation, evidence)
    r["participants"] = participants
    r["steps"] = steps
    r["outcomes"] = outcomes
    r["reaction_scale"] = F(*scale, evidence=evidence) if scale else F(evidence=evidence)
    if solvent:
        r["solvent_conditions"]["composition"] = F(solvent, solvent, evidence)
    r["procedure"]["original_text"] = F(
        f"Experimental procedure transcribed from PDF pages {', '.join(map(str, pages))}; ordered steps and charges are represented explicitly.",
        f"Experimental procedure on PDF pages {', '.join(map(str, pages))}.", evidence)
    r["procedure"]["addition_order"] = F(
        [s["operation"]["normalized_value"] for s in steps],
        [s["step_id"] for s in steps], evidence)
    r["notes"] = [F(n, n, evidence) for n in notes]
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved)}
    r["completeness"] = "partial"
    if classification is not None:
        r["classification"] = classification
    return r


def add_simple(mid, name, pages, smiles=None, chirality="achiral", geometry="not_applicable",
               representation="encoded", note=None):
    evidence = ev(*pages)
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{"document_id": DOC, "compound_label": F(name, name, evidence),
                         "evidence_ids": evidence}]
    m["names"]["original"] = F(name, name, evidence)
    m["names"]["english"] = F(name, name, evidence)
    m["names"]["abbreviations"] = F(evidence=evidence)
    m["names"]["aliases"] = F(evidence=evidence)
    m["cas_numbers"] = F(evidence=evidence)
    if smiles:
        diag = st.stereo_inspect(smiles)
        iso = diag["canonical_isomeric_smiles"]
        m["representation_type"] = "smiles"
        m["smiles"] = F({"source_name": name, "candidate": smiles}, iso, evidence,
                         method="structure_derivation", verification="inferred",
                         derivation="Name-derived machine candidate; parsing does not establish source identity.")
        m["stereochemistry"]["isomeric_smiles"] = F(smiles, iso, evidence,
                         method="structure_derivation", verification="inferred",
                         derivation="RDKit canonical isomeric SMILES for the name-derived candidate.")
        m["stereochemistry"]["identity_key"] = F(None, diag["stereo_graph_key"], evidence,
                         method="structure_derivation", verification="inferred",
                         derivation="Stereo graph key for the machine candidate, not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F(evidence=evidence)
        m["generic_template"] = None
        m["stereochemistry"]["isomeric_smiles"] = F(evidence=evidence)
        m["stereochemistry"]["identity_key"] = F(evidence=evidence)
    m["stereochemistry"]["atom_centers"] = F(evidence=evidence)
    m["stereochemistry"]["double_bonds"] = F(evidence=evidence)
    m["stereochemistry"]["axial"] = F(evidence=evidence)
    m["stereochemistry"]["mixture"] = F(evidence=evidence)
    m["stereochemistry"]["assignment"] = F(evidence=evidence)
    for key, val in (("chirality_status", chirality), ("geometry_status", geometry),
                     ("representation_status", representation if smiles else "not_encoded")):
        m["stereochemistry"][key] = F(val, val, evidence, method="inference",
             verification="inferred", derivation="Status assigned from the explicit name/structure type; machine record remains needs_review.")
    text = note or ("Name-derived candidate structure; source identity remains needs_review." if smiles
                    else "Source-named material retained without an unsafe complete SMILES.")
    m["identity_note"] = F(text, text, evidence, method="inference", verification="inferred", derivation=text)
    m["external_properties"] = []
    return m


def refine_dense(mid, exact_name, pages, centers=None, double=None, mixture=None,
                 chirality="specified", geometry="not_applicable"):
    m = deepcopy(OLD_M[mid])
    evidence = list(dict.fromkeys(m["occurrences"][0]["evidence_ids"] + ev(*pages)))
    label = m["occurrences"][0]["compound_label"]["normalized_value"] or mid
    m["occurrences"] = [{"document_id": DOC, "compound_label": F(label, label, evidence),
                         "evidence_ids": evidence}]
    m["names"]["original"] = F(exact_name, exact_name, ev(*pages))
    m["names"]["english"] = F(exact_name, exact_name, ev(*pages))
    if centers is not None:
        m["stereochemistry"]["atom_centers"] = F(centers, centers, ev(*pages))
    if double is not None:
        m["stereochemistry"]["double_bonds"] = F(double, double, ev(*pages))
    if mixture is not None:
        m["stereochemistry"]["mixture"] = F(mixture, mixture, ev(*pages))
    m["stereochemistry"]["chirality_status"] = F(chirality, chirality, ev(*pages),
        method="inference", verification="inferred", derivation="Status follows explicit configurational locants in the source heading/drawing.")
    m["stereochemistry"]["geometry_status"] = F(geometry, geometry, ev(*pages),
        method="inference", verification="inferred", derivation="Status follows explicit E/Z locants or the source drawing; no missing geometry is treated as achiral.")
    m["stereochemistry"]["representation_status"] = F("not_encoded", "not_encoded", ev(*pages),
        method="inference", verification="inferred", derivation="Dense source structure and its stereochemistry were not converted to a complete isomeric SMILES.")
    m["identity_note"] = F(
        "Exact source heading and stereochemical locants retained; complete SMILES remains empty because the dense drawing was not independently encoded.",
        "Exact source name retained; dense structure remains unencoded and needs_review.", ev(*pages),
        method="inference", verification="inferred",
        derivation="Avoids asserting a chemically complete graph from name parsing alone.")
    return m


new_molecules = {}
new_molecules["NH4Cl"] = add_simple("NH4Cl", "ammonium chloride", (7, 8, 9), "[NH4+].[Cl-]")
new_molecules["nBuLi"] = add_simple("nBuLi", "n-butyllithium", (7, 9), None,
                                      note="Organolithium formulation retained by source name; aggregation/solvation and a complete representation are unresolved.")
new_molecules["N2"] = add_simple("N2", "nitrogen", (7,), "N#N")
new_molecules["NaHCO3"] = add_simple("NaHCO3", "sodium bicarbonate", (8, 9, 10), "[Na+].O=C([O-])O")
new_molecules["Na2S2O3"] = add_simple("Na2S2O3", "sodium thiosulfate", (8, 9), None,
                                       note="Aqueous sodium thiosulfate is source-named; concentration and exact solution composition are unreported.")
new_molecules["silica_gel"] = add_simple("silica_gel", "silica gel", (7, 8, 9, 10), None,
                                          chirality="not_applicable", geometry="not_applicable",
                                          note="Chromatography stationary phase; variable material, not represented as a unique molecule.")
new_molecules["minus_methoxy_ipc_borane"] = add_simple(
    "minus_methoxy_ipc_borane", "(−)-methoxydiisopinocampheylborane", (7, 9), None,
    chirality="specified", geometry="not_applicable", representation="not_encoded",
    note="Chiral reagent configuration is source-specified by its (−) name; no complete stereochemical graph is asserted.")
new_molecules["minus_methoxy_ipc_borane"]["stereochemistry"]["assignment"] = F(
    "(−)-methoxydiisopinocampheylborane", {"source_assignment": "(−) reagent", "independent_RS_assignment": False},
    ev(7, 9), method="original_extraction")
new_molecules["phosphate_adduct_pair"] = add_simple(
    "phosphate_adduct_pair", "inseparable diastereomeric alpha-phosphate Michael-adduct pair", (11,), None,
    chirality="mixture", geometry="specified", representation="not_encoded",
    note="Tentatively assigned source-drawn pair of diastereomers; full graph and individual configurations are not independently encoded.")
pm = new_molecules["phosphate_adduct_pair"]
pm["stereochemistry"]["atom_centers"] = F(
    {"source": "footnote-39 drawing retains the precursor stereotetrad; newly generated diastereomeric relationship is not assigned"},
    {"specified_precursor_centers": True, "new_center_assignment": "unknown"}, ev(11))
pm["stereochemistry"]["double_bonds"] = F(
    {"source": "footnote-39 drawing shows the enone and terminal alkene"},
    {"retained_descriptively": True}, ev(11))
pm["stereochemistry"]["mixture"] = F(
    "a pair of inseparable diastereomers", {"type": "diastereomeric_pair", "ratio": None}, ev(11))
pm["stereochemistry"]["assignment"] = F(
    "tentatively assigned as the Michael adduct of alpha-phosphate to the acrylate based on 1H NMR and LS-ESI",
    {"assignment": "tentative alpha-phosphate Michael adduct", "basis": ["1H NMR", "LS-ESI"]}, ev(11))


dense_specs = {
    "cmp_18": ("(2Z,4E,7R,8R,9R,10S)-7-((tert-butyldimethylsilyl)oxy)-11-((tert-butyldiphenylsilyl)oxy)-4-ethyl-9-methoxy-8,10-dimethylundeca-2,4-dienoic acid", (7, 11), ["7R", "8R", "9R", "10S"], ["2Z", "4E"], "specified"),
    "cmp_21": ("(5R,6R)-5-ethyl-6-(((2R,3S,4S,5R)-4-methoxy-3,5-dimethyl-6-oxotetrahydro-2H-pyran-2-yl)methyl)-5,6-dihydro-2H-pyran-2-one", (7,), ["5R", "6R", "side-chain 2R,3S,4S,5R"], None, "not_applicable"),
    "cmp_22": ("(5S,6R,7R,8S)-5-allyl-7-methoxy-2,2,3,3,6,8,12,12-octamethyl-11,11-diphenyl-4,10-dioxa-3,11-disilatridecane", (7,), ["5S", "6R", "7R", "8S"], None, "not_applicable"),
    "cmp_23": ("(3S,4R,5R,6S)-3-((tert-butyldimethylsilyl)oxy)-7-((tert-butyldiphenylsilyl)oxy)-5-methoxy-4,6-dimethylheptanal", (7,), ["3S", "4R", "5R", "6S"], None, "not_applicable"),
    "cmp_24": ("(3R,4R,6S,7R,8R,9S)-6-((tert-butyldimethylsilyl)oxy)-10-((tert-butyldiphenylsilyl)oxy)-3-ethyl-8-methoxy-7,9-dimethyldec-1-en-4-yl acrylate", (7, 8), ["3R", "4R", "6S", "7R", "8R", "9S"], None, "not_applicable"),
    "cmp_25": ("(5R,6R)-6-((2S,3R,4R,5S)-2-((tert-butyldimethylsilyl)oxy)-6-((tert-butyldiphenylsilyl)oxy)-4-methoxy-3,5-dimethylhexyl)-5-ethyl-5,6-dihydro-2H-pyran-2-one", (8,), ["5R", "6R", "side-chain 2S,3R,4R,5S"], None, "not_applicable"),
    "cmp_26a": ("(5R,6R)-6-((2S,3R,4R,5S)-2-((tert-butyldimethylsilyl)oxy)-6-hydroxy-4-methoxy-3,5-dimethylhexyl)-5-ethyl-5,6-dihydro-2H-pyran-2-one", (8,), ["5R", "6R", "side-chain 2S,3R,4R,5S"], None, "not_applicable"),
    "cmp_26b": ("(2Z,4E,7S,8R,9R,10S)-7-((tert-butyldimethylsilyl)oxy)-11-((tert-butyldiphenylsilyl)oxy)-4-ethyl-9-methoxy-8,10-dimethylundeca-2,4-dienoic acid", (8,), ["7S", "8R", "9R", "10S"], ["2Z", "4E"], "specified"),
    "cmp_26c": ("(1R,5R,7S,9R)-9-ethyl-7-((2S,3R,4S)-5-hydroxy-3-methoxy-4-methylpentan-2-yl)-2,6-dioxabicyclo[3.3.1]nonan-3-one", (8,), ["1R", "5R", "7S", "9R", "side-chain 2S,3R,4S"], None, "not_applicable"),
    "cmp_32": ("(5R,6R)-6-((2S,3R,4R,5S)-(E)-2-((tert-butyldimethylsilyl)oxy)-4-methoxy-3,5-dimethyl-8-oxonon-6-en-1-yl)-5-ethyl-5,6-dihydro-2H-pyran-2-one", (8, 9), ["5R", "6R", "side-chain 2S,3R,4R,5S"], ["6E"], "specified"),
    "cmp_35": ("(3R,4R,6R,7R,8R,9S)-6-((tert-butyldimethylsilyl)oxy)-10-((tert-butyldiphenylsilyl)oxy)-3-ethyl-8-methoxy-7,9-dimethyldec-1-en-4-yl cinnamate", (9,), ["3R", "4R", "6R", "7R", "8R", "9S"], {"source": "drawn cinnamate geometry retained descriptively"}, "specified"),
    "cmp_36": ("(3R,4R,6R,7R,8R,9S)-6-((tert-butyldimethylsilyl)oxy)-3-ethyl-10-hydroxy-8-methoxy-7,9-dimethyldec-1-en-4-yl cinnamate", (9,), ["3R", "4R", "6R", "7R", "8R", "9S"], {"source": "drawn cinnamate geometry retained descriptively"}, "specified"),
    "cmp_37": ("(3R,4R,6R,7R,8R,9S)-(E)-6-((tert-butyldimethylsilyl)oxy)-3-ethyl-8-methoxy-7,9-dimethyl-12-oxotrideca-1,10-dien-4-yl cinnamate", (9,), ["3R", "4R", "6R", "7R", "8R", "9S"], ["10E", "drawn cinnamate geometry"], "specified"),
    "cmp_30": ("(5R,6R)-6-((2R,3R,4R,5S)-(E)-2-((tert-butyldimethylsilyl)oxy)-4-methoxy-3,5-dimethyl-8-oxonon-6-en-1-yl)-5-ethyl-5,6-dihydro-2H-pyran-2-one", (9, 10), ["5R", "6R", "side-chain 2R,3R,4R,5S"], ["6E"], "specified"),
    "cmp_1": ("(5R,6R)-5-ethyl-6-((2R,3S,4R,5S)-(E)-2-hydroxy-4-methoxy-3,5-dimethylnon-7-en-1-yl)-5,6-dihydro-2H-pyran-2-one ((−)-pironetin)", (10,), ["5R", "6R", "side-chain 2R,3S,4R,5S"], ["7E"], "specified"),
}
for mid, (name, pages, centers, dbl, geom) in dense_specs.items():
    new_molecules[mid] = refine_dense(mid, name, pages, centers=centers, double=dbl,
                                      chirality="specified", geometry=geom)


reactions = []

# 17 -> 18 (page 7 full procedure)
e = ev(3, 4, 7)
parts = [P("substrate", "cmp_17", "reactant", e, mass=Q(25.2, "mg"), amount=Q(0.038, "mmol")),
         P("thf", "THF", "solvent", e, volume=Q(1.5, "mL")),
         P("tbaf", "TBAF", "reagent", e, volume=Q(77, "uL", "77 μL")),
         P("nh4cl", "NH4Cl", "workup_agent", e, stage="workup", volume=Q(2, "mL", "2 mL aqueous solution")),
         P("etoac", "EtOAc", "workup_agent", e, stage="workup", volume=Q(30, "mL", "3 × 10 mL")),
         P("brine", "brine", "workup_agent", e, stage="workup"),
         P("dry", "Na2SO4", "workup_agent", e, stage="workup"),
         P("silica", "silica_gel", "purification_agent", e, stage="purification"),
         P("pe", "petroleum_ether", "purification_agent", e, stage="purification"),
         P("ea", "EtOAc", "purification_agent", e, stage="purification"),
         P("product", "cmp_18", "product", e)]
steps = [S("tbaf_treatment", e, ("substrate", "thf", "tbaf"), "add 1.0 M TBAF in THF and stir", "room temperature", Q(2, "h")),
         S("quench_extract", e, ("nh4cl", "etoac", "brine", "dry"), "quench with aqueous NH4Cl; extract three times with EtOAc; wash with brine, dry, filter, and concentrate", stage="workup"),
         S("chromatography", e, ("silica", "pe", "ea", "product"), "flash chromatography on silica gel, PE/EA 15:1", stage="purification")]
reactions.append(rxn("scheme5_17_to_18", (3, 4, 7), parts, steps,
    [O("yield", "yield", ("cmp_18",), e, Q(75, "%", "19 mg, 75% yield"), "isolated after flash chromatography", "isolated")],
    "TBAF-promoted pyrone fragmentation to dienoic acid 18",
    unresolved=("The source gives TBAF volume and solution concentration but not an independently stated mmol charge.",),
    notes=("Page 7 provides the full procedure for the Scheme 5 transformation.",),
    location="Scheme 5 and page-7 compound 18 procedure",
    scale=("lactone 17, 25.2 mg (0.038 mmol)", {"number": 0.038, "unit": "mmol", "reference_molecule_id": "cmp_17"}), solvent="anhydrous THF"))

# 17 -> 20 first stage of the page-7 sequence.
e = ev(3, 4, 7)
parts = [P("substrate", "cmp_17", "reactant", e, mass=Q(5, "mg"), amount=Q(0.007, "mmol")),
         P("meoh", "MeOH", "solvent", e, volume=Q(1, "mL")), P("tsoh", "TsOH", "reagent", e, mass=Q(10, "mg")),
         P("silica", "silica_gel", "purification_agent", e, stage="purification"), P("product", "cmp_20", "product", e)]
steps = [S("acidic_desilylation", e, ("substrate", "meoh", "tsoh"), "add TsOH in one portion and heat until TLC indicated full consumption", Q(60, "°C")),
         S("short_plug", e, ("silica", "product"), "remove solvent and pass residue through a short plug of silica gel", stage="purification")]
reactions.append(rxn("scheme5_17_to_20", (3, 4, 7), parts, steps,
    [O("yield", "yield", ("cmp_20",), e, Q(80, "%", "80% (Scheme 5)"), None, None,
       "Scheme 5 step yield from 17 to 20")], "acidic removal of the two silyl groups",
    unresolved=("Page 7 calls 20 the resulting alcohol but gives no isolated mass or yield there; the 80% value is only printed on Scheme 5.",
                "Numeric reaction time and the yield measurement basis are unreported."),
    notes=("The isolated/plug-treated alcohol was immediately used in the following oxidation; its amount charged to that stage is not printed.",),
    location="Scheme 5 and first stage of page-7 compound 21 procedure",
    scale=("pyrone 17, 5 mg (0.007 mmol)", {"number": 0.007, "unit": "mmol", "reference_molecule_id": "cmp_17"}), solvent="MeOH"))

# 20 -> 21 second stage.
e = ev(3, 4, 7)
parts = [P("substrate", "cmp_20", "reactant", e), P("dcm", "CH2Cl2", "solvent", e, volume=Q(0.5, "mL")),
         P("tempo", "TEMPO", "catalyst", e, mass=Q(2, "mg"), amount=Q(0.01, "mmol")),
         P("baib", "BAIB", "reagent", e, mass=Q(6, "mg"), amount=Q(0.02, "mmol")),
         P("silica", "silica_gel", "purification_agent", e, stage="purification"),
         P("pe", "petroleum_ether", "purification_agent", e, stage="purification"),
         P("ea", "EtOAc", "purification_agent", e, stage="purification"), P("product", "cmp_21", "product", e)]
steps = [S("oxidation_lactonization", e, ("substrate", "dcm", "tempo", "baib"), "redissolve the resulting alcohol in anhydrous CH2Cl2; add TEMPO and (diacetoxyiodo)benzene respectively; stir", "room temperature", Q(2, "h", "additional 2 h")),
         S("chromatography", e, ("silica", "pe", "ea", "product"), "remove solvent and flash chromatograph on silica gel, PE/EA 1:2", stage="purification")]
reactions.append(rxn("scheme5_20_to_21", (3, 4, 7), parts, steps,
    [O("yield", "yield", ("cmp_21",), e, Q(79, "%", "1.8 mg, 79% yield"), "isolated after flash chromatography", "isolated")],
    "TEMPO/BAIB oxidation and lactonization to bislactone 21",
    unresolved=("The amount of intermediate 20 entering the oxidation is not reported, so the exact denominator for the printed 79% cannot be independently reconstructed.",),
    notes=("The page-7 prose is continuous from 17 through plug-treated alcohol 20 to 21; Scheme 5 separately depicts the 20-to-21 arrow.",),
    location="Scheme 5 and second stage of page-7 compound 21 procedure", solvent="anhydrous CH2Cl2"))

# 3a -> 22.
e = ev(4, 7)
parts = [P("substrate", "cmp_3a", "reactant", e, mass=Q(384, "mg"), amount=Q(0.931, "mmol")),
         P("cui", "CuI", "reagent", e, mass=Q(35, "mg"), amount=Q(0.186, "mmol")),
         P("vinyl", "vinylMgBr", "reagent", e, volume=Q(8.3, "mL", "8.3 mL of 0.56 M solution")),
         P("thf_cui", "THF", "solvent", e, volume=Q(8.0, "mL")), P("thf_sub", "THF", "solvent", e, volume=Q(6.0, "mL")),
         P("lutidine", "lutidine", "reagent", e, volume=Q(1.3, "mL"), amount=Q(11.2, "mmol")),
         P("tbsotf", "TBSOTf", "reagent", e, volume=Q(1.8, "mL"), amount=Q(7.45, "mmol")),
         P("nh4cl", "NH4Cl", "workup_agent", e, stage="workup", volume=Q(15, "mL", "15 mL aqueous solution")),
         P("etoac", "EtOAc", "workup_agent", e, stage="workup", volume=Q(30, "mL", "3 × 10 mL")),
         P("brine", "brine", "workup_agent", e, stage="workup"), P("dry", "Na2SO4", "workup_agent", e, stage="workup"),
         P("silica", "silica_gel", "purification_agent", e, stage="purification"), P("pe", "petroleum_ether", "purification_agent", e, stage="purification"),
         P("ea", "EtOAc", "purification_agent", e, stage="purification"), P("product", "cmp_22", "product", e)]
steps = [S("vinyl_reagent_addition", e, ("cui", "thf_cui", "vinyl"), "add vinylmagnesium bromide solution dropwise to CuI in anhydrous THF", Q(-20, "°C")),
         S("pre_stir", e, ("cui", "vinyl", "thf_cui"), "stir at the addition temperature", Q(-20, "°C"), Q(30, "min")),
         S("substrate_addition", e, ("substrate", "thf_sub"), "add the epoxide 3a solution"),
         S("warm_to_zero", e, ("substrate", "vinyl", "cui"), "allow the mixture to warm to 0 °C; monitor by TLC", Q(0, "°C")),
         S("silylation_addition", e, ("lutidine", "tbsotf"), "after completion add 2,6-lutidine and TBSOTf"),
         S("silylation_hold", e, ("substrate", "lutidine", "tbsotf"), "stir until the alcohol intermediate was completely consumed", Q(0, "°C")),
         S("workup", e, ("nh4cl", "etoac", "brine", "dry"), "slowly quench with aqueous NH4Cl; separate, extract, wash with brine, dry, filter, and concentrate", stage="workup"),
         S("chromatography", e, ("silica", "pe", "ea", "product"), "flash chromatography on silica gel, PE/EA 100:1", stage="purification")]
reactions.append(rxn("scheme6_3a_to_22", (4, 7), parts, steps,
    [O("yield", "yield", ("cmp_22",), e, Q(95, "%", "490 mg, 95% yield"), "isolated after flash chromatography", "isolated")],
    "vinylative epoxide opening and in-pot TBS protection of the C7 epimer",
    unresolved=("The duration of warming to 0 °C and the duration of the TBS-protection hold are not numerically reported.",),
    location="Scheme 6 and page-7 compound 22 procedure",
    scale=("epoxide 3a, 384 mg (0.931 mmol)", {"number": 0.931, "unit": "mmol", "reference_molecule_id": "cmp_3a"}), solvent="anhydrous THF"))

# 22 -> 23.
e = ev(3, 4, 7)
parts = [P("substrate", "cmp_22", "reactant", e, mass=Q(408, "mg"), amount=Q(0.736, "mmol")),
         P("dcm", "CH2Cl2", "solvent", e, volume=Q(16, "mL")), P("ozone", "ozone", "reagent", e),
         P("nitrogen", "N2", "reagent", e), P("pph3", "PPh3", "reagent", e, mass=Q(0.97, "g"), amount=Q(3.68, "mmol")),
         P("silica", "silica_gel", "purification_agent", e, stage="purification"), P("pe", "petroleum_ether", "purification_agent", e, stage="purification"),
         P("ea", "EtOAc", "purification_agent", e, stage="purification"), P("product", "cmp_23", "product", e)]
steps = [S("cool", e, ("substrate", "dcm"), "cool the alkene solution", Q(-78, "°C")),
         S("ozonolysis", e, ("substrate", "ozone", "dcm"), "treat with ozone until the solution turned light blue", Q(-78, "°C")),
         S("nitrogen_flush", e, ("nitrogen",), "flush with nitrogen gas", Q(-78, "°C"), Q(10, "min")),
         S("reductive_workup", e, ("pph3",), "add PPh3 in one portion and stir", "room temperature", Q(2, "h", "additional 2 h")),
         S("chromatography", e, ("silica", "pe", "ea", "product"), "remove solvent and flash chromatograph on silica gel, PE/EA 60:1", stage="purification")]
reactions.append(rxn("scheme6_22_to_23", (3, 4, 7), parts, steps,
    [O("yield", "yield", ("cmp_23",), e, Q(87, "%", "356 mg, 87% yield"), "isolated after flash chromatography", "isolated")],
    "ozonolysis/PPh3 workup to aldehyde 23", location="Scheme 6 and page-7 compound 23 procedure",
    scale=("alkene 22, 408 mg (0.736 mmol)", {"number": 0.736, "unit": "mmol", "reference_molecule_id": "cmp_22"}), solvent="CH2Cl2"))


def allylation_acylation_record(suffix, pages, substrate_mid, product_mid, product_name,
                               pentene, borane, bf3, substrate, naoh, peroxide,
                               acyl_mid, acyl_charge, base_charge, dmap=None, acyl_temp="room temperature",
                               acyl_time=None, yield_value=None, conflict=False):
    evidence = ev(*pages)
    parts = [P("substrate", substrate_mid, "reactant", evidence, mass=Q(*substrate[:2]), amount=Q(*substrate[2:])),
             P("pentene", "cis_2_pentene", "reagent", evidence, volume=Q(*pentene[:2]), amount=Q(*pentene[2:])),
             P("tbuok", "tBuOK", "reagent", evidence, volume=Q(1.4, "mL", "1.4 mL of 1.0 M in THF")),
             P("thf_initial", "THF", "solvent", evidence, volume=Q(1.5, "mL")),
             P("nbuli", "nBuLi", "reagent", evidence, volume=Q(0.59, "mL", "0.59 mL of 2.4 M in hexanes")),
             P("borane", "minus_methoxy_ipc_borane", "reagent", evidence, mass=Q(*borane[:2]), amount=Q(*borane[2:])),
             P("thf_borane", "THF", "solvent", evidence, volume=Q(3.5, "mL")),
             P("bf3", "BF3_Et2O", "reagent", evidence, volume=Q(*bf3[:2]), amount=Q(*bf3[2:])),
             P("naoh", "NaOH", "reagent", evidence, volume=Q(*naoh[:2])),
             P("peroxide", "H2O2", "reagent", evidence, volume=Q(*peroxide[:2])),
             P("etoac_extract", "EtOAc", "workup_agent", evidence, stage="workup", volume=Q(30, "mL", "3 × 10 mL")),
             P("brine", "brine", "workup_agent", evidence, stage="workup"), P("dry", "Na2SO4", "workup_agent", evidence, stage="workup"),
             P("dcm_acyl", "CH2Cl2", "solvent", evidence, volume=Q(30 if suffix.endswith("23_to_24") else 20, "mL")),
             P("base", "Et3N", "reagent", evidence, volume=Q(*base_charge[:2]), amount=Q(*base_charge[2:])),
             P("acyl", acyl_mid, "reagent", evidence, mass=Q(*acyl_charge[:2]) if acyl_charge[1] == "mg" else None,
               volume=Q(*acyl_charge[:2]) if acyl_charge[1] != "mg" else None, amount=Q(*acyl_charge[2:])),
             P("silica", "silica_gel", "purification_agent", evidence, stage="purification"),
             P("pe", "petroleum_ether", "purification_agent", evidence, stage="purification"), P("ea", "EtOAc", "purification_agent", evidence, stage="purification"),
             P("product", product_mid, "product", evidence)]
    if dmap:
        parts.insert(-4, P("dmap", "DMAP", "reagent", evidence, mass=Q(*dmap[:2]), amount=Q(*dmap[2:])))
    steps = [S("nbuli_addition", evidence, ("pentene", "tbuok", "thf_initial", "nbuli"), "slowly add nBuLi to cis-2-pentene/potassium tert-butoxide suspension", Q(-78, "°C")),
             S("minus50_hold", evidence, ("pentene", "tbuok", "nbuli"), "warm and stir", Q(-50, "°C"), Q(10, "min")),
             S("recool", evidence, ("pentene", "tbuok", "nbuli"), "re-cool", Q(-78, "°C")),
             S("borane_addition", evidence, ("borane", "thf_borane"), "add the chiral methoxyborane solution", Q(-78, "°C")),
             S("borane_hold", evidence, ("borane",), "stir", Q(-78, "°C"), Q(30, "min")),
             S("bf3_substrate_addition", evidence, ("bf3", "substrate"), "add BF3 etherate and then the aldehyde, respectively", Q(-78, "°C")),
             S("allylation_hold", evidence, ("substrate", "bf3", "borane"), "stir until starting material was consumed by TLC", Q(-78, "°C"), Q(6, "h")),
             S("oxidative_workup", evidence, ("naoh", "peroxide"), "slowly add aqueous NaOH and 30 wt% H2O2, then heat at reflux", "reflux", Q(1, "h"), stage="workup"),
             S("extract", evidence, ("etoac_extract", "brine", "dry"), "extract three times with EtOAc, wash with brine, dry, and concentrate; use residue directly", stage="workup")]
    acyl_pids = ["base", "acyl", "dcm_acyl"] + (["dmap"] if dmap else [])
    steps += [S("acylation_addition", evidence, tuple(acyl_pids), f"add base{', DMAP' if dmap else ''} and {product_name} acylating reagent", Q(0, "°C")),
              S("acylation_hold", evidence, tuple(acyl_pids + ["product"]), "warm and stir", acyl_temp, acyl_time),
              S("chromatography", evidence, ("silica", "pe", "ea", "product"),
                "remove solvent and flash chromatograph on silica gel, " + ("PE/EA 15:1" if suffix.endswith("23_to_24") else "PE/EA 150:1"), stage="purification")]
    outs = []
    if conflict:
        o = st.new_outcome("yield", "yield", (product_mid,), evidence)
        o["value"] = F(["68% (experimental procedure and Scheme 10)", "78% (page-5 prose)"], None, evidence,
                       status="conflict", alternatives=[
                           F("361 mg, 68% yield", {"number": 68, "unit": "%", "raw_text": "361 mg, 68% yield"}, ev(9)),
                           F("78%", {"number": 78, "unit": "%", "raw_text": "78%"}, ev(5))])
        o["measurement_basis"] = F("isolated after flash chromatography for the 68% procedure value; basis of prose 78% not reported", None, evidence, status="conflict",
                                    alternatives=[F("afforded 35 (361 mg, 68% yield) after flash chromatography", "isolated after flash chromatography", ev(9)), F("obtained in good overall yield (78%)", None, ev(5))])
        o["yield_type"] = F("isolated/overall yield conflict", "conflict", evidence)
        outs.append(o)
    else:
        outs.append(O("yield", "yield", (product_mid,), evidence, yield_value, "isolated after flash chromatography", "isolated"))
    unresolved = ["Ramp durations for warming to −50 °C, re-cooling to −78 °C, and the final warm-up are unreported."]
    notes = []
    if conflict:
        unresolved.append("Compound 35 yield conflicts: 68% in Scheme 10 and the detailed page-9 procedure versus 78% in page-5 prose; neither value is selected.")
        notes.append("Both source values and their separate evidence are retained as a conflict.")
    return rxn(suffix, pages, parts, steps, outs,
               "stereoselective pentenylation, oxidative workup, and " + ("acrylation" if acyl_mid == "acryloyl_chloride" else "cinnamoylation"),
               unresolved=tuple(unresolved), notes=tuple(notes),
               location=("Scheme 6 and pages 7-8 compound 24 procedure" if suffix.endswith("23_to_24") else "Scheme 10, page-5 prose, and page-9 compound 35 procedure"),
               scale=(f"{substrate_mid}, {substrate[2]} {substrate[3]}", {"number": substrate[2], "unit": substrate[3], "reference_molecule_id": substrate_mid}),
               solvent="THF for allylation; CH2Cl2 for acylation")


reactions.append(allylation_acylation_record(
    "scheme6_23_to_24", (3, 4, 7, 8), "cmp_23", "cmp_24", "acryloyl",
    (0.46, "mL", 4.27, "mmol"), (562, "mg", 1.78, "mmol"), (0.22, "mL", 1.78, "mmol"),
    (396, "mg", 0.71, "mmol"), (1.2, "mL"), (1.2, "mL"),
    "acryloyl_chloride", (0.45, "mL", 5.5, "mmol"), (1.0, "mL", 7.4, "mmol"),
    acyl_temp="room temperature", acyl_time=Q(3, "h", "another 3 h"),
    yield_value=Q(65, "%", "315 mg, 65% yield")))

# 24 -> 25.
e = ev(4, 8)
parts = [P("substrate", "cmp_24", "reactant", e, mass=Q(137, "mg"), amount=Q(0.20, "mmol")),
         P("dcm_sub", "CH2Cl2", "solvent", e, volume=Q(40, "mL")),
         P("grubbs", "GrubbsII", "catalyst", e, mass=Q(34, "mg"), amount=Q(0.04, "mmol")),
         P("dcm_cat", "CH2Cl2", "solvent", e, volume=Q(3, "mL")),
         P("silica", "silica_gel", "purification_agent", e, stage="purification"), P("pe", "petroleum_ether", "purification_agent", e, stage="purification"),
         P("ea", "EtOAc", "purification_agent", e, stage="purification"), P("product", "cmp_25", "product", e)]
steps = [S("catalyst_addition", e, ("substrate", "dcm_sub", "grubbs", "dcm_cat"), "add Grubbs-II solution to the ester solution; monitor by TLC", Q(40, "°C")),
         S("chromatography", e, ("silica", "pe", "ea", "product"), "remove solvent and flash chromatograph on silica gel, PE/EA 25:1", stage="purification")]
reactions.append(rxn("scheme6_24_to_25", (4, 8), parts, steps,
    [O("isolated_yield", "yield", ("cmp_25",), e, Q(80, "%", "105 mg, 80% yield"), "isolated after flash chromatography", "isolated"),
     O("brsm_yield", "yield", ("cmp_25",), e, Q(91, "%", "brsm: 91%"), "based on recovered starting material; recovered amount not reported", "yield based on recovered starting material")],
    "ring-closing metathesis to pyrone 25",
    unresolved=("Reaction time is reported only by TLC completion; the recovered starting-material amount underlying 91% brsm is unreported.",),
    location="Scheme 6 and page-8 compound 25 procedure",
    scale=("ester 24, 137 mg (0.20 mmol)", {"number": 0.20, "unit": "mmol", "reference_molecule_id": "cmp_24"}), solvent="freshly distilled CH2Cl2"))

# 25 -> 26a/26b/26c.
e = ev(4, 8)
parts = [P("substrate", "cmp_25", "reactant", e, mass=Q(56, "mg"), amount=Q(0.085, "mmol")),
         P("thf", "THF", "solvent", e, volume=Q(4, "mL")), P("tbaf", "TBAF", "reagent", e, volume=Q(0.17, "mL", "0.17 mL of 1.0 M in THF")),
         P("nh4cl", "NH4Cl", "workup_agent", e, stage="workup", volume=Q(1.0, "mL", "1.0 mL aqueous solution")),
         P("etoac", "EtOAc", "workup_agent", e, stage="workup", volume=Q(30, "mL", "3 × 10 mL")), P("brine", "brine", "workup_agent", e, stage="workup"),
         P("dry", "Na2SO4", "workup_agent", e, stage="workup"), P("silica", "silica_gel", "purification_agent", e, stage="purification"),
         P("pe", "petroleum_ether", "purification_agent", e, stage="purification"), P("ea", "EtOAc", "purification_agent", e, stage="purification"),
         P("p26a", "cmp_26a", "product", e), P("p26b", "cmp_26b", "product", e), P("p26c", "cmp_26c", "product", e)]
steps = [S("tbaf_addition", e, ("substrate", "thf", "tbaf"), "add TBAF solution dropwise", "ambient temperature"),
         S("tbaf_hold", e, ("substrate", "thf", "tbaf"), "stir", "ambient temperature", Q(2, "h")),
         S("workup", e, ("nh4cl", "etoac", "brine", "dry"), "add aqueous NH4Cl; separate, extract three times, wash, dry, filter, and concentrate", stage="workup"),
         S("chromatography", e, ("silica", "pe", "ea", "p26a", "p26b", "p26c"), "flash chromatography on silica gel, PE/EA 5:1", stage="purification")]
reactions.append(rxn("scheme6_25_desilylation_branches", (4, 8), parts, steps,
    [O("yield_26a", "yield", ("cmp_26a",), e, Q(37, "%", "13 mg, 37% yield"), "isolated after flash chromatography", "isolated"),
     O("yield_26b", "yield", ("cmp_26b",), e, Q(20, "%", "11 mg, 20% yield"), "isolated after flash chromatography", "isolated"),
     O("yield_26c", "yield", ("cmp_26c",), e, Q(32, "%", "8.3 mg, 32% yield"), "isolated after flash chromatography", "isolated")],
    "TBAF desilylation with fragmentation and bridged oxa-Michael branches",
    unresolved=("The source does not report the fate of the remaining mass balance.",),
    location="Scheme 6 and page-8 three-product TBAF procedure",
    scale=("pyrone 25, 56 mg (0.085 mmol)", {"number": 0.085, "unit": "mmol", "reference_molecule_id": "cmp_25"}), solvent="freshly distilled THF"))


def oxidation_hwe(suffix, pages, substrate_mid, product_mid, substrate_mass, substrate_mmol,
                  dcm_vol, phosphonate_vol, phosphonate_mmol, thf_base, base_vol,
                  thf_crude, hold_time, yield_mass, yield_pct, temp_conflict=False):
    evidence = ev(*pages)
    parts = [P("substrate", substrate_mid, "reactant", evidence, mass=Q(substrate_mass, "mg"), amount=Q(substrate_mmol, "mmol")),
             P("dcm", "CH2Cl2", "solvent", evidence, volume=Q(dcm_vol, "mL")),
             P("dmp", "DMP", "reagent", evidence, mass=Q(60, "mg"), amount=Q(0.14, "mmol")),
             P("nahco3", "NaHCO3", "workup_agent", evidence, stage="workup"), P("thiosulfate", "Na2S2O3", "workup_agent", evidence, stage="workup"),
             P("etoac", "EtOAc", "workup_agent", evidence, stage="workup", volume=Q(30, "mL", "3 × 10 mL")),
             P("brine", "brine", "workup_agent", evidence, stage="workup"), P("dry", "Na2SO4", "workup_agent", evidence, stage="workup"),
             P("phosphonate", "phosphonate", "reagent", evidence, volume=Q(phosphonate_vol, "uL", f"{phosphonate_vol} μL"), amount=Q(phosphonate_mmol, "mmol")),
             P("thf_base", "THF", "solvent", evidence, volume=Q(thf_base, "mL")),
             P("tbuok", "tBuOK", "reagent", evidence, volume=Q(base_vol, "uL", f"{base_vol} μL of 1.0 M in THF")),
             P("thf_crude", "THF", "solvent", evidence, volume=Q(thf_crude, "mL")),
             P("water", "water", "workup_agent", evidence, stage="workup", volume=("a few drops", "a few drops")),
             P("silica", "silica_gel", "purification_agent", evidence, stage="purification"), P("pe", "petroleum_ether", "purification_agent", evidence, stage="purification"),
             P("ea", "EtOAc", "purification_agent", evidence, stage="purification"), P("product", product_mid, "product", evidence)]
    steps = [S("dmp_addition", evidence, ("substrate", "dcm", "dmp"), "add Dess-Martin periodinane in one portion", Q(0, "°C")),
             S("post_addition_hold", evidence, ("substrate", "dmp"), "hold after addition", Q(0, "°C"), Q(10, "min")),
             S("oxidation_hold", evidence, ("substrate", "dmp"), "warm to room temperature and stir", "room temperature", Q(2, "h")),
             S("oxidation_workup", evidence, ("nahco3", "thiosulfate", "etoac", "brine", "dry"), "quench with saturated NaHCO3 and Na2S2O3 solutions; extract, wash, dry, and concentrate; use crude directly", stage="workup"),
             S("base_addition", evidence, ("phosphonate", "thf_base", "tbuok"), "slowly add 1.0 M potassium tert-butoxide to dimethyl acetylmethylphosphonate solution", Q(0, "°C")),
             S("anion_hold", evidence, ("phosphonate", "tbuok"), "stir", "room temperature", Q(1, "h")),
             S("crude_addition", evidence, ("thf_crude",), "add the crude aldehyde solution dropwise", Q(0, "°C"))]
    if temp_conflict:
        temp_fact = F(["room temperature (page-9 detailed procedure)", "45 °C (page-5 prose)"], None, evidence, status="conflict",
                      alternatives=[F("room temperature", "room temperature", ev(9)), F("45 °C", {"number": 45, "unit": "°C", "raw_text": "45 °C"}, ev(5))])
        hwe = S("hwe_hold", evidence, ("phosphonate", "tbuok", "substrate", "product"), "warm and stir overnight", time=hold_time)
        hwe["temperature"] = temp_fact
        steps.append(hwe)
    else:
        steps.append(S("hwe_hold", evidence, ("phosphonate", "tbuok", "substrate", "product"), "warm to room temperature and stir", "room temperature", hold_time))
    steps += [S("water_quench", evidence, ("water",), "quench with a few drops of water", stage="workup"),
              S("chromatography", evidence, ("silica", "pe", "ea", "product"), "remove solvent and purify on silica gel, " + ("PE/EA 5:1" if product_mid == "cmp_32" else "PE/EA 10:1"), stage="purification")]
    unresolved = ["The crude aldehyde amount is not independently measured before the HWE stage."]
    notes = []
    if temp_conflict:
        unresolved.append("HWE temperature conflicts: page-5 prose says 45 °C, while the detailed page-9 procedure says room temperature overnight.")
        notes.append("Both temperature statements are retained as evidence-linked alternatives without selecting one.")
    return rxn(suffix, pages, parts, steps,
               [O("yield", "yield", (product_mid,), evidence, Q(yield_pct, "%", f"{yield_mass} mg, {yield_pct}%"), "isolated after silica-gel purification", "isolated")],
               "Dess-Martin oxidation followed by Horner-Wadsworth-Emmons homologation",
               unresolved=tuple(unresolved), notes=tuple(notes),
               location=("Scheme 9 and pages 8-9 compound 32 procedure" if product_mid == "cmp_32" else "Scheme 10, page-5 prose, and page-9 compound 37 procedure"),
               scale=(f"{substrate_mid}, {substrate_mass} mg ({substrate_mmol} mmol)", {"number": substrate_mmol, "unit": "mmol", "reference_molecule_id": substrate_mid}),
               solvent="CH2Cl2 for oxidation; anhydrous THF for HWE")


reactions.append(oxidation_hwe("scheme9_26a_to_32", (4, 5, 8, 9), "cmp_26a", "cmp_32", 14.7, 0.035, 5, 36, 0.26, 3.0, 215, 2.0, Q(4, "h"), 13.3, 84))

reactions.append(allylation_acylation_record(
    "scheme10_14_to_35", (5, 9), "cmp_14", "cmp_35", "cinnamoyl",
    (0.45, "mL", 4.23, "mmol"), (558, "mg", 1.76, "mmol"), (217, "uL", 1.76, "mmol"),
    (393, "mg", 0.71, "mmol"), (1.1, "mL"), (1.1, "mL"),
    "cinnamoyl_chloride", (847, "mg", 5.1, "mmol"), (1.1, "mL", 7.6, "mmol"),
    dmap=(310, "mg", 2.54, "mmol"), acyl_temp=Q(40, "°C"), acyl_time=Q(5, "h", "another 5 h"), conflict=True))

# 35 -> 36.
e = ev(5, 9)
parts = [P("substrate", "cmp_35", "reactant", e, mass=Q(46.6, "mg"), amount=Q(0.06, "mmol")), P("thf", "THF", "solvent", e, volume=Q(3.0, "mL")),
         P("tbaf", "TBAF", "reagent", e, volume=Q(185, "uL", "185 μL of 1.0 M in THF")),
         P("nh4cl", "NH4Cl", "workup_agent", e, stage="workup", volume=("a few drops of aqueous solution", "a few drops")),
         P("etoac", "EtOAc", "workup_agent", e, stage="workup", volume=Q(30, "mL", "3 × 10 mL")), P("brine", "brine", "workup_agent", e, stage="workup"),
         P("dry", "Na2SO4", "workup_agent", e, stage="workup"), P("silica", "silica_gel", "purification_agent", e, stage="purification"),
         P("pe", "petroleum_ether", "purification_agent", e, stage="purification"), P("ea", "EtOAc", "purification_agent", e, stage="purification"), P("product", "cmp_36", "product", e)]
steps = [S("tbaf_addition", e, ("substrate", "thf", "tbaf"), "slowly add TBAF solution", "room temperature"),
         S("overnight_hold", e, ("substrate", "tbaf", "thf"), "stir overnight", "room temperature", ("overnight", "overnight")),
         S("workup", e, ("nh4cl", "etoac", "brine", "dry"), "quench with a few drops aqueous NH4Cl; extract, wash, dry, filter, and concentrate", stage="workup"),
         S("chromatography", e, ("silica", "pe", "ea", "product"), "flash chromatography on silica gel, PE/EA 15:1", stage="purification")]
reactions.append(rxn("scheme10_35_to_36", (5, 9), parts, steps,
    [O("yield", "yield", ("cmp_36",), e, Q(82, "%", "25.7 mg, 82%"), "isolated after flash chromatography", "isolated")],
    "selective TBDPS deprotection to primary alcohol 36", location="Scheme 10 and page-9 compound 36 procedure",
    scale=("ester 35, 46.6 mg (0.06 mmol)", {"number": 0.06, "unit": "mmol", "reference_molecule_id": "cmp_35"}), solvent="freshly distilled THF"))

reactions.append(oxidation_hwe("scheme10_36_to_37", (5, 9), "cmp_36", "cmp_37", 19.4, 0.037, 3, 38, 0.27, 2.0, 227, 2.0, ("overnight", "overnight"), 15.8, 78, temp_conflict=True))

# 37 -> 30.
e = ev(5, 9)
parts = [P("substrate", "cmp_37", "reactant", e, mass=Q(15.8, "mg"), amount=Q(0.028, "mmol")), P("toluene_sub", "toluene", "solvent", e, volume=Q(5, "mL")),
         P("grubbs", "GrubbsII", "catalyst", e, mass=Q(4.8, "mg"), amount=Q(0.0057, "mmol")), P("toluene_cat", "toluene", "solvent", e, volume=Q(0.5, "mL")),
         P("silica", "silica_gel", "purification_agent", e, stage="purification"), P("pe", "petroleum_ether", "purification_agent", e, stage="purification"),
         P("ea", "EtOAc", "purification_agent", e, stage="purification"), P("product", "cmp_30", "product", e)]
steps = [S("heat", e, ("substrate", "toluene_sub"), "warm the enone solution", Q(100, "°C")),
         S("catalyst_addition", e, ("grubbs", "toluene_cat"), "slowly add Grubbs-II solution; monitor by TLC", Q(100, "°C")),
         S("chromatography", e, ("silica", "pe", "ea", "product"), "remove solvent and flash chromatograph on silica gel, PE/EA 5:1", stage="purification")]
reactions.append(rxn("scheme10_37_to_30", (5, 9), parts, steps,
    [O("yield", "yield", ("cmp_30",), e, Q(80, "%", "10.1 mg, 80% yield"), "isolated after flash chromatography", "isolated")],
    "ring-closing metathesis to enone-pyrone 30",
    unresolved=("Reaction duration is reported only by TLC monitoring.",),
    location="Scheme 10 and page-9 compound 30 procedure",
    scale=("enone 37, 15.8 mg (0.028 mmol)", {"number": 0.028, "unit": "mmol", "reference_molecule_id": "cmp_37"}), solvent="freshly distilled toluene"))

# 30 -> pironetin 1.
e = ev(4, 5, 9, 10)
parts = [P("substrate", "cmp_30", "reactant", e, mass=Q(8.0, "mg"), amount=Q(0.017, "mmol")), P("meoh", "MeOH", "solvent", e, volume=Q(2.0, "mL")),
         P("tosylhydrazide", "TsNHNH2", "reagent", e, mass=Q(6.6, "mg"), amount=Q(0.035, "mmol")),
         P("borohydride", "NaBH4", "reagent", e, mass=Q(14, "mg"), amount=Q(0.40, "mmol")), P("acoh", "AcOH", "solvent", e, volume=Q(2, "mL", "2 mL glacial acetic acid")),
         P("hcl", "HCl", "reagent", e), P("nahco3", "NaHCO3", "workup_agent", e, stage="workup", volume=Q(5, "mL", "5 mL of 2 N aqueous solution")),
         P("etoac", "EtOAc", "workup_agent", e, stage="workup", volume=Q(30, "mL", "3 × 10 mL")), P("brine", "brine", "workup_agent", e, stage="workup"),
         P("dry", "Na2SO4", "workup_agent", e, stage="workup"), P("silica", "silica_gel", "purification_agent", e, stage="purification"),
         P("pe", "petroleum_ether", "purification_agent", e, stage="purification"), P("ea", "EtOAc", "purification_agent", e, stage="purification"), P("product", "cmp_1", "product", e)]
steps = [S("hydrazone_formation", e, ("substrate", "meoh", "tosylhydrazide"), "add 4-methylbenzenesulfonhydrazide in one portion and warm until complete consumption", Q(60, "°C")),
         S("solvent_removal", e, ("substrate", "tosylhydrazide"), "remove methanol under reduced pressure"),
         S("borohydride_addition", e, ("borohydride", "acoh"), "add sodium borohydride solution in glacial acetic acid dropwise"),
         S("reduction_hold", e, ("borohydride", "acoh"), "stir; monitor by TLC", "room temperature", Q(3, "h")),
         S("acidification", e, ("hcl",), "add concentrated aqueous HCl to adjust pH to 2-3"),
         S("acidic_desilylation", e, ("hcl", "substrate"), "stir the acidified solution", Q(70, "°C"), Q(3, "h")),
         S("neutralization", e, ("nahco3",), "slowly add 2 N aqueous NaHCO3", Q(0, "°C"), stage="workup"),
         S("extraction", e, ("etoac", "brine", "dry"), "extract three times with EtOAc; wash with brine, dry, and concentrate", stage="workup"),
         S("chromatography", e, ("silica", "pe", "ea", "product"), "flash chromatography on silica gel, PE/EA 3:1", stage="purification")]
reactions.append(rxn("scheme8_30_to_pironetin_onepot", (4, 5, 9, 10), parts, steps,
    [O("overall_yield", "yield", ("cmp_1",), e, Q(76, "%", "4.2 mg, 76% yield"), "isolated after flash chromatography", "isolated overall yield",
       "overall from pyrone 30 through hydrazone/reduction/rearrangement and acidic desilylation to (−)-pironetin 1")],
    "hydrazone formation, NaBH4 reduction/diazene rearrangement, and acidic TBS removal",
    unresolved=("The hydrazone-formation duration is reported only as complete consumption, not a numeric time.",
                "The printed 76% is a start-to-end multistage yield and is not assigned to the final desilylation alone."),
    notes=("The final product is reported as white solid and identical with literature data; characterization is source context, not a separate synthesis record.",),
    location="Scheme 8/Scheme 10 context and pages 9-10 pironetin procedure",
    scale=("pyrone 30, 8.0 mg (0.017 mmol)", {"number": 0.017, "unit": "mmol", "reference_molecule_id": "cmp_30"}), solvent="MeOH followed by glacial acetic acid"))

# Footnote 39 belongs to the existing 28 -> 29 HWE experiment; augment rather than duplicate.
base = deepcopy(OLD_R[f"{DOC}:scheme8_28_to_29"])
e = ev(4, 5, 11)
base["evidence_ids"] = e
base["location"] = F("Scheme 8, adjacent discussion, and footnote 39", {"pages": [4, 5, 11], "region": "28 to 29 HWE and phosphate Michael-adduct footnote"}, e)
base["participants"].append(P("phosphate_adduct_pair", "phosphate_adduct_pair", "product", ev(11)))
base["outcomes"].append(O("phosphate_adduct_yield", "yield", ("phosphate_adduct_pair",), ev(11),
                           Q(64, "%", "64% yield"), "isolated; tentative structure based on 1H NMR and LS-ESI", "isolated"))
base["notes"] = [F("Footnote 39 reports an inseparable pair of diastereomers tentatively assigned as the alpha-phosphate Michael adduct to the acrylate.",
                    "Footnote-39 phosphate Michael-adduct pair retained as an observed side product.", ev(11))]
base["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                  "unresolved": ["The detailed preparation of 28 and 29 is assigned to prior work ref. 16 and is not reproduced on pages 7-11.",
                                 "The side-product pair is tentative and inseparable; individual configurations and ratio are unreported.",
                                 "The 47% single-pass, 80% after two cycles, 86% brsm scheme value, and 64% side-product result have different scopes and must not be collapsed."]}
# Add an observed C-C branch in addition to the retained DB-C-C HWE label.
changes = list(base["classification"]["bond_changes"])
cid = f"b{len(changes)+1}"
changes.append(F(None, {"change_id": cid, "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
                        "before_order": 0, "after_order": 1, "change_type": "formed", "branch": "footnote-39 phosphate Michael adduct"},
                 ev(11), method="structure_derivation", verification="inferred",
                 derivation="The footnote-39 drawing shows a new C-C bond from the phosphonate alpha carbon to the acrylate beta carbon; atom mapping is not asserted."))
base["classification"]["bond_changes"] = changes
labels = deepcopy(base["classification"]["labels"]["normalized_value"] or [])
if not any(x["code"] == "C-C" for x in labels):
    labels.append({"category_number": 2, "code": "C-C", "supporting_change_ids": [cid]})
base["classification"]["labels"] = F(None, labels, e, method="structure_derivation", verification="inferred",
    derivation="Multi-label assignment covers the main HWE C=C formation and the observed footnote-39 Michael-addition C-C branch.")
base["classification"]["coverage"] = F("partial", "partial", e, method="inference", verification="inferred",
    derivation="Dense structures are not atom mapped; only source-visible net bond formations are indexed.")
reactions.append(base)


# Assemble only the delta plus every molecule referenced by that delta.
required_mids = set()
for r in reactions:
    required_mids.update(p["molecule_id"] for p in r["participants"])
molecules = []
for mid in sorted(required_mids):
    molecules.append(new_molecules.get(mid) or deepcopy(OLD_M[mid]))

bundle = st.empty(DATASET["documents"])
bundle["molecules"] = molecules
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": "pages7_10_characterization_blocks", "document_id": DOC, "pages": [7, 8, 9, 10],
     "evidence_ids": ev(7, 8, 9, 10), "kind": "characterization",
     "summary": "Optical rotations, 1H/13C NMR, FT-IR and HRMS blocks following compounds 18, 21-26, 32, 35-37, 30 and pironetin 1 are characterization, not independent synthesis experiments."},
    {"context_id": "page10_back_matter_and_unlinked_si", "document_id": DOC, "pages": [10],
     "evidence_ids": ev(10), "kind": "source_boundary",
     "summary": "Competing-interest declaration, acknowledgments, and an online Supporting Information link. No SI file was supplied or associated; the link is not treated as linked SI evidence."},
    {"context_id": "pages10_11_references", "document_id": DOC, "pages": [10, 11],
     "evidence_ids": ev(10, 11), "kind": "references",
     "summary": "Reference list and cited prior literature are not current-paper experiments. Author footnotes 31 and 39 are separately retained where they add current-paper stereochemical or side-product information."},
    {"context_id": "page11_footnote31_compound18_geometry", "document_id": DOC, "pages": [11],
     "evidence_ids": ev(11), "kind": "stereochemical_note",
     "summary": "Footnote 31 states that the C2-C3 stereochemistry of dienoic acid 18 remains Z and refers to unavailable SI for details; Z is retained, but no SI evidence is claimed."},
]

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 7, "state": "attempted", "regions": [
            {"region": "left column, compound 18 and compound 21 procedures", "finding": "full TBAF fragmentation and two-stage desilylation/oxidation procedures; spectral blocks separated as characterization", "evidence_ids": ev(7)},
            {"region": "lower left through right column, compound 22 and compound 23 procedures", "finding": "full vinylative epoxide opening/TBS protection and ozonolysis procedures", "evidence_ids": ev(7)},
            {"region": "lower right column, beginning of compound 24 procedure", "finding": "low-temperature pentenylation setup continues on page 8", "evidence_ids": ev(7)}],
         "record_ids": [f"{DOC}:scheme5_17_to_18", f"{DOC}:scheme5_17_to_20", f"{DOC}:scheme5_20_to_21", f"{DOC}:scheme6_3a_to_22", f"{DOC}:scheme6_22_to_23", f"{DOC}:scheme6_23_to_24"],
         "unresolved": ["Several TLC-defined endpoints have no numeric duration; page-7 characterization is not converted to synthesis records."]},
        {"page": 8, "state": "attempted", "regions": [
            {"region": "upper left, continuation and completion of compound 24 procedure", "finding": "oxidative workup, acrylation, purification and 65% isolated yield; merged with page-7 start", "evidence_ids": ev(8)},
            {"region": "middle left, compound 25 procedure", "finding": "RCM experiment with 80% isolated and 91% brsm outcomes", "evidence_ids": ev(8)},
            {"region": "lower left and upper right, TBAF treatment of 25 and 26a/26b/26c characterization", "finding": "one three-product experiment with separate 37%, 20%, and 32% isolated yields", "evidence_ids": ev(8)},
            {"region": "lower right, beginning of compound 32 procedure", "finding": "Dess-Martin oxidation and HWE setup continues on page 9", "evidence_ids": ev(8)}],
         "record_ids": [f"{DOC}:scheme6_23_to_24", f"{DOC}:scheme6_24_to_25", f"{DOC}:scheme6_25_desilylation_branches", f"{DOC}:scheme9_26a_to_32"],
         "unresolved": ["The 91% brsm denominator lacks a reported recovered-starting-material amount; spectral blocks remain source context."]},
        {"page": 9, "state": "attempted", "regions": [
            {"region": "upper left, completion/characterization of compound 32", "finding": "84% isolated enone 32; merged with page-8 procedure", "evidence_ids": ev(9)},
            {"region": "left column, compounds 35 and 36 procedures", "finding": "full cinnamate sequence and selective deprotection; compound 35 retains 68% versus 78% source conflict", "evidence_ids": ev(9)},
            {"region": "right column, compounds 37 and 30 procedures", "finding": "oxidation/HWE and RCM experiments; HWE temperature conflict retained between page-5 prose and page-9 procedure", "evidence_ids": ev(9)},
            {"region": "bottom/right continuation into page 10", "finding": "compound 30 characterization and pironetin procedure begin/continue across page boundary", "evidence_ids": ev(9)}],
         "record_ids": [f"{DOC}:scheme9_26a_to_32", f"{DOC}:scheme10_14_to_35", f"{DOC}:scheme10_35_to_36", f"{DOC}:scheme10_36_to_37", f"{DOC}:scheme10_37_to_30", f"{DOC}:scheme8_30_to_pironetin_onepot"],
         "unresolved": ["Page-5 prose says 45 °C for the HWE to 37, but the detailed page-9 procedure says room temperature overnight; neither is selected."]},
        {"page": 10, "state": "attempted", "regions": [
            {"region": "upper left, completion of pironetin 1 procedure and characterization", "finding": "multistage preparation gives 4.2 mg, 76% isolated (−)-pironetin; spectra/rotation are characterization", "evidence_ids": ev(10)},
            {"region": "lower left, declarations/acknowledgments/Supplementary data link", "finding": "non-experimental back matter; SI link noted but no SI associated", "evidence_ids": ev(10)},
            {"region": "right column, references 1-19", "finding": "cited prior literature, not current-paper experiments", "evidence_ids": ev(10)}],
         "record_ids": [f"{DOC}:scheme8_30_to_pironetin_onepot"],
         "unresolved": ["Hydrazone formation ends at complete consumption without a numeric time; online SI was not supplied and is not associated."]},
        {"page": 11, "state": "attempted", "regions": [
            {"region": "reference list and footnote 31", "finding": "prior literature plus current-paper note that compound 18 retains Z C2-C3 geometry; cited SI details unavailable", "evidence_ids": ev(11)},
            {"region": "footnote 39 and displayed structure", "finding": "real observed 64% isolated inseparable diastereomeric alpha-phosphate Michael-adduct pair; merged into the existing 28-to-29 HWE experiment", "evidence_ids": ev(11)},
            {"region": "remaining references 40-43", "finding": "cited prior work, not extracted as current-paper experiments", "evidence_ids": ev(11)}],
         "record_ids": [f"{DOC}:scheme5_17_to_18", f"{DOC}:scheme8_28_to_29"],
         "unresolved": ["Footnote-39 side-product assignment is tentative; individual diastereomer configurations and ratio are unreported; footnote-31 SI evidence is unavailable."]},
    ]
}

raw_inventory = {
    "document_id": DOC,
    "pages_reviewed_from_original_images": [7, 8, 9, 10, 11],
    "updated_reaction_ids": [r["reaction_id"] for r in reactions],
    "new_molecule_ids": sorted(set(new_molecules) - set(OLD_M)),
    "source_context_ids": [x["context_id"] for x in bundle["source_context"]],
    "chemical_accuracy_verified": False,
    "notes": ["Same experiments shown earlier as scheme arrows were enriched with page-7-to-10 procedures rather than duplicated.",
              "Page-11 footnote 39 was merged into the existing HWE experiment as a side-product branch.",
              "All records remain needs_review; no SI association was inferred."]
}

st.write(RUN / "submitted.json", bundle)
st.write(RUN / "raw_model/pages_0007_0011_extraction.json", raw_inventory)
st.write(RUN / "coverage/pages_0007_0011.json", coverage)
print(json.dumps({"reactions": len(reactions), "molecules_in_bundle": len(molecules),
                  "new_molecules": len(raw_inventory["new_molecule_ids"]), "pages": 5}, ensure_ascii=False))
