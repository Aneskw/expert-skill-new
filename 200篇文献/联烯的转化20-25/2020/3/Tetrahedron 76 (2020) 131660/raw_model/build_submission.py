#!/usr/bin/env python3
"""Build v3.1 machine records for Gu et al., Tetrahedron 76 (2020) 131660, PDF pp. 1-6."""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import standard_records as st

DOC = "doc_c724145241bb10a8"
RUN = ROOT / "runs/current"
DATASET = st.read(RUN / "dataset.json")
EV = {p: [f"{DOC}:page{p}:image", f"{DOC}:page{p}:text"] for p in range(1, 7)}


def ev(*pages):
    return list(dict.fromkeys(x for p in pages for x in EV[p]))


def F(raw=None, norm=None, evidence=(), method="original_extraction", status=None,
      verification=None, derivation=None, alternatives=()):
    return st.fact(raw, norm, evidence, method=method, status=status,
                   verification=verification, derivation=derivation,
                   alternatives=alternatives)


def inferred(value, evidence, why):
    return F(value, value, evidence, method="inference", verification="inferred", derivation=why)


bundle = st.empty(DATASET["documents"])
molecules = {}
stereo_checks = {}


def molecule(mid, label, pages, name=None, smiles=None, ch="unresolved", geom="unresolved",
             rep="unresolved", centers=None, double=None, axial=None, assignment=None, note=None):
    if mid in molecules:
        return mid
    evidence = ev(*pages)
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{"document_id": DOC, "compound_label": F(label, label, evidence), "evidence_ids": evidence}]
    source_name = name or label
    m["names"]["original"] = F(source_name, source_name, evidence)
    m["names"]["english"] = F(source_name, source_name, evidence)
    m["names"]["abbreviations"] = F(evidence=evidence)
    m["names"]["aliases"] = F(evidence=evidence)
    m["cas_numbers"] = F(evidence=evidence)
    if smiles:
        d = st.stereo_inspect(smiles)
        stereo_checks[mid] = d
        canonical = d["canonical_isomeric_smiles"]
        m["representation_type"] = "smiles"
        m["smiles"] = F({"source_name": source_name, "candidate": smiles}, canonical, evidence,
            method="structure_derivation", verification="inferred",
            derivation="Candidate graph derived from the explicit source name; parsing does not verify source identity.")
        m["stereochemistry"]["isomeric_smiles"] = F({"candidate": smiles}, canonical, evidence,
            method="structure_derivation", verification="inferred",
            derivation="RDKit canonical isomeric SMILES for the machine candidate graph.")
        m["stereochemistry"]["identity_key"] = F(None, d["stereo_graph_key"], evidence,
            method="structure_derivation", verification="inferred",
            derivation="RDKit stereo graph key; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F(evidence=evidence)
        m["generic_template"] = None
        m["stereochemistry"]["isomeric_smiles"] = F(evidence=evidence)
        m["stereochemistry"]["identity_key"] = F(evidence=evidence)
    m["stereochemistry"]["atom_centers"] = F(centers, centers, evidence) if centers is not None else F(evidence=evidence)
    m["stereochemistry"]["double_bonds"] = F(double, double, evidence) if double is not None else F(evidence=evidence)
    m["stereochemistry"]["axial"] = F(axial, axial, evidence) if axial is not None else F(evidence=evidence)
    m["stereochemistry"]["assignment"] = F(assignment, assignment, evidence) if assignment is not None else F(evidence=evidence)
    m["stereochemistry"]["mixture"] = F(evidence=evidence)
    m["stereochemistry"]["chirality_status"] = inferred(ch, evidence, "Assigned from source name/drawing; missing labels were not interpreted as achiral or racemic.")
    m["stereochemistry"]["geometry_status"] = inferred(geom, evidence, "Assigned from source name/drawing and candidate graph where applicable.")
    m["stereochemistry"]["representation_status"] = inferred(rep, evidence, "Reports whether source-visible stereochemistry is encoded in the stored representation.")
    identity_note = note or ("Name-derived machine candidate; source identity remains needs_review." if smiles else "Dense source drawing retained by label; no unsafe complete SMILES asserted.")
    m["identity_note"] = F(identity_note, identity_note, evidence, method="inference", verification="inferred", derivation=identity_note)
    m["external_properties"] = []
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


# Dense paper-labelled structures. Relative/absolute wedges and named E/Z information
# are retained descriptively; complete isomeric SMILES are intentionally not guessed.
dense = {
    "cmp_1": ("1", (1, 2, 5), "pironetin", "specified", "specified"),
    "cmp_2": ("2", (2,), "stereotetrad 2", "specified", "unresolved"),
    "cmp_3": ("3", (2, 3), "epoxide 3", "specified", "not_applicable"),
    "cmp_3a": ("3a", (2, 4), "C7-epimer 3a", "specified", "not_applicable"),
    "cmp_4": ("4", (2, 3, 6), "syn-2,3-allenol 4", "specified", "unspecified"),
    "cmp_5": ("5", (2, 3, 6), "syn,anti,syn-triol 5", "specified", "not_applicable"),
    "cmp_6": ("6", (2,), "(S)-Roche ester 6", "specified", "not_applicable"),
    "cmp_7": ("7", (2,), "TBDPS-protected Roche ester 7", "specified", "not_applicable"),
    "cmp_8": ("8", (2,), "Brown propargylborane 8", "partially_specified", "not_applicable"),
    "cmp_9": ("9", (2, 3, 6), "desired epoxide 9", "specified", "not_applicable"),
    "cmp_10": ("10", (2, 3), "C7-epimer epoxide 10", "specified", "not_applicable"),
    "alcohol_9": ("alcohol derived from major epoxide 9", (2, 6), "alcohol derived from major epoxide 9", "specified", "unspecified"),
    "cmp_R11": ("(R)-MTPA-11", (2, 6), "(R)-MTPA ester 11", "specified", "unspecified"),
    "cmp_S11": ("(S)-MTPA-11", (6,), "(S)-MTPA ester 11", "specified", "unspecified"),
    "cmp_12": ("12", (3,), "beta-hydroxy aldehyde 12", "specified", "not_applicable"),
    "cmp_13": ("13", (3,), "alkene 13", "specified", "unspecified"),
    "cmp_14": ("14", (3, 5), "aldehyde 14", "specified", "not_applicable"),
    "cmp_15": ("15", (3, 5), "Z-pentenyl diisopinocampheylborane 15", "partially_specified", "specified"),
    "cmp_16": ("16", (3, 4, 5, 6), "acrylate 16", "specified", "unspecified"),
    "cmp_17": ("17", (3, 4, 6), "pyrone 17", "specified", "specified"),
    "cmp_18": ("18", (3, 4), "dienoic acid 18", "specified", "specified"),
    "cmp_19": ("19", (3, 4), "mono-desilylated alcohol 19 (not observed)", "specified", "specified"),
    "cmp_20": ("20", (3, 4), "diol 20", "specified", "specified"),
    "cmp_20a": ("20a", (3, 4), "lactol 20a", "specified", "specified"),
    "cmp_21": ("21", (3, 4), "bislactone 21", "specified", "specified"),
    "cmp_22": ("22", (4,), "alkene 22", "specified", "unspecified"),
    "cmp_23": ("23", (3, 4), "aldehyde 23", "specified", "not_applicable"),
    "cmp_24": ("24", (4,), "acrylate 24", "specified", "unspecified"),
    "cmp_25": ("25", (3, 4), "pyrone 25", "specified", "specified"),
    "cmp_26a": ("26a", (4, 5), "alcohol 26a", "specified", "specified"),
    "cmp_26b": ("26b", (4,), "dienic acid 26b", "specified", "specified"),
    "cmp_26c": ("26c", (4,), "caged compound 26c", "specified", "not_applicable"),
    "cmp_27": ("27", (4,), "desired iodination product 27 (not observed)", "specified", "specified"),
    "cmp_27a": ("27a", (4,), "tosylate 27a", "specified", "specified"),
    "cmp_28": ("28", (4, 5), "primary alcohol 28", "specified", "unspecified"),
    "cmp_29": ("29", (4, 5), "E-enone 29", "specified", "specified"),
    "cmp_30": ("30", (4, 5), "enone-pyrone 30", "specified", "specified"),
    "cmp_31": ("31", (4, 5), "TBS-protected pironetin 31", "specified", "specified"),
    "cmp_32": ("32", (4, 5), "7-epi E-enone 32", "specified", "specified"),
    "cmp_33": ("33", (4, 5), "TBS-protected 7-epi framework 33", "specified", "specified"),
    "cmp_34": ("34", (5,), "double-bond-migrated crotonyl product 34", "specified", "partially_specified"),
    "cmp_35": ("35", (5,), "cinnamate 35", "specified", "specified"),
    "cmp_36": ("36", (5,), "primary alcohol 36", "specified", "specified"),
    "cmp_37": ("37", (5,), "E-enone cinnamate 37", "specified", "specified"),
    "cmp_7epi": ("7-epi-pironetin", (5,), "7-epi-pironetin (not obtained)", "specified", "specified"),
    "side_thf": ("five-membered side product", (3,), "five-membered cyclized side product", "unresolved", "unresolved"),
}
for mid, (label, pages, name, ch, geom) in dense.items():
    centers = {"source": "wedge/dash or configurational name", "converted_to_independent_RS": False} if ch in {"specified", "partially_specified"} else None
    db = {"source": "drawn alkene geometry or E/Z name", "retained_descriptively": True} if geom in {"specified", "partially_specified"} else None
    axial = {"allene": "terminal allene; no R_a/S_a assignment reported"} if mid == "cmp_4" else None
    molecule(mid, label, pages, name=name, ch=ch, geom=geom, rep="not_encoded", centers=centers, double=db, axial=axial,
             assignment={"source_assignment": name, "independent_assignment": False} if centers or db else None)


simple = {
    "THF": ("tetrahydrofuran", "C1CCOC1"), "CH2Cl2": ("dichloromethane", "ClCCl"),
    "MeOH": ("methanol", "CO"), "toluene": ("toluene", "Cc1ccccc1"), "hexane": ("hexane", "CCCCCC"),
    "acetone": ("acetone", "CC(C)=O"),
    "Et3N": ("triethylamine", "CCN(CC)CC"), "DMAP": ("4-dimethylaminopyridine", "CN(C)c1ccncc1"),
    "DMF": ("N,N-dimethylformamide", "CN(C)C=O"), "oxalyl_chloride": ("oxalyl chloride", "O=C(Cl)C(=O)Cl"),
    "Et2O": ("diethyl ether", "CCOCC"), "water": ("water", "O"), "NaOH": ("sodium hydroxide", "[Na+].[OH-]"),
    "H2O2": ("hydrogen peroxide", "OO"), "EtOAc": ("ethyl acetate", "CCOC(C)=O"),
    "Na2SO4": ("sodium sulfate", "[Na+].[Na+].[O-]S(=O)(=O)[O-]"), "MeI": ("methyl iodide", "CI"),
    "NaH": ("sodium hydride", "[Na+].[H-]"), "TsCl": ("p-toluenesulfonyl chloride", "Cc1ccc(S(=O)(=O)Cl)cc1"),
    "NsCl": ("o-nitrobenzenesulfonyl chloride", "O=[N+]([O-])c1ccccc1S(=O)(=O)Cl"),
    "anisyl_SO2Cl": ("4-methoxybenzenesulfonyl chloride", "COc1ccc(S(=O)(=O)Cl)cc1"),
    "DEAD": ("diethyl azodicarboxylate", "CCOC(=O)N=NC(=O)OCC"),
    "DIAD": ("diisopropyl azodicarboxylate", "CC(C)OC(=O)N=NC(=O)OC(C)C"),
    "PPh3": ("triphenylphosphine", "P(c1ccccc1)(c1ccccc1)c1ccccc1"), "CuI": ("copper(I) iodide", "[Cu]I"),
    "lutidine": ("2,6-lutidine", "Cc1cccc(C)n1"), "ozone": ("ozone", "[O-][O+]=O"),
    "acryloyl_chloride": ("acryloyl chloride", "C=CC(=O)Cl"), "TsOH": ("p-toluenesulfonic acid", "Cc1ccc(S(=O)(=O)O)cc1"),
    "TEMPO": ("TEMPO", "CC1(C)CCCC(C)(C)N1[O]"), "I2": ("iodine", "II"), "NaI": ("sodium iodide", "[Na+].[I-]"),
    "phosphonate": ("dimethyl acetylmethylphosphonate", "CC(=O)CP(=O)(OC)OC"),
    "tBuOK": ("potassium tert-butoxide", "[K+].CC(C)(C)[O-]"),
    "TsNHNH2": ("p-toluenesulfonhydrazide", "Cc1ccc(S(=O)(=O)NN)cc1"),
    "NaBH4": ("sodium borohydride", "[Na+].[BH4-]"), "AcOH": ("acetic acid", "CC(=O)O"),
    "HCl": ("hydrochloric acid", "Cl"), "EtOH": ("ethanol", "CCO"),
    "cis_2_pentene": ("cis-2-pentene", "CC/C=C\\C"),
    "cinnamoyl_chloride": ("cinnamoyl chloride", "O=C(Cl)/C=C/c1ccccc1"),
    "two_methyl_two_butene": ("2-methyl-2-butene", "CC=C(C)C"),
}
for mid, (name, smi) in simple.items():
    d = st.stereo_inspect(smi)
    atom_unspec = any(x["type"].startswith("Atom_") and x["specified"] == "Unspecified" for x in d["potential_stereo"])
    any_spec = any(x["specified"] == "Specified" for x in d["potential_stereo"])
    any_unspec = any(x["specified"] == "Unspecified" for x in d["potential_stereo"])
    molecule(mid, name, (2, 3, 4, 5, 6), name=name, smiles=smi,
             ch="unspecified" if atom_unspec else "achiral",
             geom="specified" if any_spec else "unspecified" if any_unspec else "not_applicable",
             rep="partial" if any_unspec else "encoded",
             double={"source":"name-derived E/Z or cis geometry","toolkit":d["double_bonds"]} if any_spec else None,
             assignment={"candidate_graph":d["canonical_isomeric_smiles"],"source_verified":False} if any_spec else None)


unresolved = [
    ("TBDPSCl", "tert-butyldiphenylsilyl chloride", (2,), "achiral", "not_applicable"),
    ("DIBALH", "diisobutylaluminum hydride", (2,), "achiral", "not_applicable"),
    ("Sia2BH", "disiamylborane (Sia2BH)", (3, 6), "achiral", "not_applicable"),
    ("IpcBH2", "(+)-diisopinocampheylborane", (3, 6), "specified", "not_applicable"),
    ("Ipc_complex", "((+)-IpcBH2)2·TMEDA", (6,), "specified", "not_applicable"),
    ("BH3_DMS", "borane-dimethyl sulfide complex", (6,), "achiral", "not_applicable"),
    ("BF3_Et2O", "boron trifluoride etherate", (3, 5, 6), "achiral", "not_applicable"),
    ("Trisimid", "2,4,6-triisopropylbenzenesulfonyl imidazole (Trisimid)", (2, 3), "achiral", "not_applicable"),
    ("TrisCl", "2,4,6-triisopropylbenzenesulfonyl chloride (TrisCl)", (3,), "achiral", "not_applicable"),
    ("KHMDS", "potassium hexamethyldisilazide", (3,), "achiral", "not_applicable"),
    ("vinylMgBr", "vinylmagnesium bromide", (3, 4), "achiral", "not_applicable"),
    ("TBSOTf", "tert-butyldimethylsilyl triflate", (3, 4), "achiral", "not_applicable"),
    ("lithiated_bisboronate", "lithiated di-[B(pin)2]-methane reagent", (3,), "unresolved", "not_applicable"),
    ("GrubbsII", "Grubbs second-generation catalyst", (3, 4, 5, 6), "unresolved", "not_applicable"),
    ("TBAF", "tetrabutylammonium fluoride", (3, 4, 5), "achiral", "not_applicable"),
    ("Ph3PI2", "triphenylphosphine-iodine reagent", (4,), "achiral", "not_applicable"),
    ("DMP", "Dess-Martin periodinane", (4, 5), "unresolved", "not_applicable"),
    ("crotonoyl_chloride", "crotonoyl chloride", (5,), "unspecified", "unspecified"),
    ("R_MTPA", "(R)-2-methoxy-2-phenyl-3,3,3-trifluoropropanoic acid", (6,), "specified", "not_applicable"),
    ("S_MTPA", "(S)-2-methoxy-2-phenyl-3,3,3-trifluoropropanoic acid", (6,), "specified", "not_applicable"),
    ("CH2Br2", "dibromomethane internal standard", (3,), "achiral", "not_applicable"),
    ("BAIB", "(diacetoxyiodo)benzene", (4,), "unresolved", "not_applicable"),
    ("petroleum_ether", "petroleum ether (composition not specified)", (6,), "mixture", "not_applicable"),
    ("brine", "brine (composition/concentration not specified)", (6,), "mixture", "not_applicable"),
    ("Celite", "Celite", (6,), "mixture", "not_applicable"),
]
for mid, name, pages, ch, geom in unresolved:
    molecule(mid, name, pages, name=name, ch=ch, geom=geom, rep="not_encoded",
             centers={"source":"configurational reagent name","independent_assignment":False} if ch == "specified" else None,
             assignment={"source_assignment":name,"independent_assignment":False} if ch == "specified" else None)


def quantity(number, unit, raw=None, qualifier=None):
    value = {"number": number, "unit": unit, "raw_text": raw or f"{number} {unit}"}
    if qualifier:
        value["qualifier"] = qualifier
    return raw or value["raw_text"], value


def participant(pid, mid, role, evidence, stage="reaction", **charges):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    for key, spec in charges.items():
        if key == "reference":
            continue
        if spec is None:
            continue
        raw, norm = spec
        p["charges"][key] = F(raw, norm, evidence)
    if any(charges.get(k) is not None for k in ("equivalents", "mol_percent")):
        p["equivalence_reference"] = F(charges.get("reference"), charges.get("reference"), evidence)
    return p


def step(sid, evidence, pids=(), operation=None, temperature=None, time=None,
         atmosphere=None, stage="reaction"):
    s = st.new_step(sid, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    if operation is not None:
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


def outcome(oid, kind, targets, evidence, value=None, basis=None, yield_type=None,
            scope=None, conflict_alternatives=()):
    o = st.new_outcome(oid, kind, targets, evidence)
    if conflict_alternatives:
        o["value"] = F([x[0] for x in conflict_alternatives], None, evidence, status="conflict",
                       alternatives=[F(raw, norm, evi) for raw, norm, evi in conflict_alternatives])
    elif value is not None:
        raw, norm = value if isinstance(value, tuple) else (value, value)
        o["value"] = F(raw, norm, evidence)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, evidence)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, evidence)
    if scope is not None:
        o["scope"] = F(scope, scope, evidence)
    return o


def classification(evidence, changes=(), coverage="partial", reason=None):
    c = {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": F(evidence=evidence),
        "bond_changes": [],
        "labels": None,
        "coverage": inferred(coverage, evidence, reason or "Partial bond-change assignment from source drawings; atom mapping was not asserted."),
        "limitations": ["Source drawings were visually inspected, but dense structures were not converted to atom-mapped SMILES.",
                        "Only explicitly asserted formed/promoted bonds are indexed; omitted changes remain unresolved."],
        "model_votes": [],
    }
    labels = {}
    for i, change in enumerate(changes, 1):
        a, b, before, after, kind = change
        payload = {"change_id": f"b{i}", "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                   "before_order": before, "after_order": after, "change_type": kind}
        if kind.startswith("hydrogen_"):
            payload["hydrogen_count"] = 1
        c["bond_changes"].append(F(None, payload, evidence, method="structure_derivation",
                                   verification="inferred", derivation="Coarse bond change read from the reaction drawing."))
        if kind in ("formed", "order_increased", "hydrogen_gain"):
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], dict(lab, supporting_change_ids=[]))["supporting_change_ids"].append(f"b{i}")
    c["labels"] = F(None, list(labels.values()), evidence, method="structure_derivation",
                    verification="inferred", derivation="Taxonomy labels calculated from the asserted partial bond changes.")
    return c


def reaction(suffix, pages, record_type, categories, parts, steps, outcomes=(),
             transformation=None, notes=(), unresolved=(), changes=(), location=None,
             solvent=None, scale=None, procedure=None, classification_coverage="partial"):
    evidence = ev(*pages)
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = F(", ".join(categories), list(categories), evidence) if categories else F(evidence=evidence)
    r["location"] = F(location or f"PDF page(s) {', '.join(map(str, pages))}",
                        {"pages": list(pages), "region": location or suffix}, evidence)
    if transformation:
        r["transformation_type"] = F(transformation, transformation, evidence)
    r["participants"] = parts
    r["steps"] = steps
    r["outcomes"] = list(outcomes)
    if solvent:
        r["solvent_conditions"]["composition"] = F(solvent, solvent, evidence)
    if scale:
        raw, norm = scale
        r["reaction_scale"] = F(raw, norm, evidence)
    if procedure:
        r["procedure"]["original_text"] = F(procedure, procedure, evidence)
    r["classification"] = classification(evidence, changes, classification_coverage)
    r["notes"] = [F(n, n, evidence) for n in notes]
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved)}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


def scheme_record(suffix, pages, title, parts, notes=(), changes=(), record_type="generic_scheme"):
    return reaction(suffix, pages, record_type, (), parts, [], transformation=title, notes=notes,
                    unresolved=("Scheme-level design/mechanistic context; not a standalone current-paper experiment.",),
                    changes=changes, classification_coverage="not_applicable" if not changes else "partial")


# Page 1: title/abstract/biological context only. Page 2: planning schemes and early route.
e2 = ev(2)
scheme_record("scheme1_initial_plan", (2,), "initial synthetic plan",
              [participant("p2", "cmp_2", "product", e2), participant("p3", "cmp_3", "reactant", e2),
               participant("p4", "cmp_4", "reactant", e2), participant("p5", "cmp_5", "product", e2)],
              notes=("Retrosynthetic/forward plan; arrows summarize intended access to the stereotetrad and are not independent procedures.",))

r = reaction("scheme2_6_to_7", (2,), "specific_experiment", ("other",),
    [participant("s6", "cmp_6", "reactant", e2), participant("tbdps", "TBDPSCl", "reagent", e2),
     participant("imid", "Trisimid", "reagent", e2), participant("solv", "CH2Cl2", "solvent", e2),
     participant("p7", "cmp_7", "product", e2)],
    [step("protection", e2, ("s6","tbdps","imid","solv","p7"), "TBDPS protection in CH2Cl2")],
    [outcome("yield_7", "yield", ("cmp_7",), e2, quantity(99, "%", "99%"), basis=None)],
    transformation="silyl protection", unresolved=("Yield measurement basis is not reported on the scheme.",),
    changes=(("Si","O",0,1,"formed"),), location="Scheme 2, first arrow")

r = reaction("scheme2_7_to_4", (2,), "specific_experiment", ("other",),
    [participant("s7", "cmp_7", "reactant", e2), participant("dibal", "DIBALH", "reagent", e2),
     participant("borane8", "cmp_8", "reagent", e2), participant("solv1", "CH2Cl2", "solvent", e2),
     participant("solv2", "THF", "solvent", e2), participant("p4", "cmp_4", "product", e2)],
    [step("reduction", e2, ("s7","dibal","solv1"), "add DIBAL-H", quantity(-78,"°C","−78 °C")),
     step("allenylation", e2, ("borane8","solv2","p4"), "apply in-situ Brown propargylborane 8", quantity(-100,"°C","−100 °C"))],
    [outcome("combined_yield_4", "yield", ("cmp_4",), e2, quantity(62,"%","62%"), basis=None,
             scope="combined yield over reduction and allenylation"),
     outcome("dr_4", "selectivity", ("cmp_4",), e2, ("dr 13:1", {"metric":"dr","ratio":[13,1],"order":["major syn-2,3-allenol 4","minor diastereomer"],"method":None}))],
    transformation="two-step reduction/allenylation", unresolved=("62% is a two-step combined yield, not a yield for either individual arrow.",
    "Diastereomer structures beyond the major source drawing are not assigned; method for dr is not stated here."),
    changes=(("C","C",0,1,"formed"),), location="Scheme 2, 7 to 4")

# Scheme 2 allene double hydroboration has a full page-6 procedure; merge both regions.
e26 = ev(2,6)
triol_parts = [
    participant("allene4", "cmp_4", "reactant", e26, amount=quantity(3.0,"mmol","3.0 mmol")),
    participant("sia", "Sia2BH", "reagent", e26), participant("ipc", "IpcBH2", "reagent", e26,
                equivalents=quantity(3.5,"equiv","~3.5 eq"),reference="cmp_4"),
    participant("alkene", "two_methyl_two_butene", "reagent", e26, volume=quantity(2.5,"mL","2.5 mL"), amount=quantity(23.1,"mmol","23.1 mmol")),
    participant("bh3", "BH3_DMS", "reagent", e26, volume=quantity(1.06,"mL","1.06 mL"), amount=quantity(10.5,"mmol","10.5 mmol")),
    participant("ipc_complex", "Ipc_complex", "reagent", e26, mass=quantity(2.55,"g","2.55 g"), amount=quantity(6.0,"mmol","6.0 mmol")),
    participant("bf3", "BF3_Et2O", "reagent", e26, volume=quantity(1.45,"mL","1.45 mL"), amount=quantity(11.7,"mmol","11.7 mmol")),
    participant("methanol", "MeOH", "reagent", e26, volume=quantity(0.9,"mL","0.9 mL")),
    participant("base", "NaOH", "reagent", e26, volume=quantity(12.0,"mL","12.0 mL")),
    participant("peroxide", "H2O2", "reagent", e26, volume=quantity(12.0,"mL","12.0 mL")),
    participant("thf_sia", "THF", "solvent", e26, volume=quantity(9.3,"mL","9.3 mL")),
    participant("thf_ipc", "THF", "solvent", e26, volume=quantity(12,"mL","12 mL")),
    participant("thf_allene", "THF", "solvent", e26, volume=quantity(6,"mL","6.0 mL")),
    participant("etoac_add", "EtOAc", "workup_agent", e26, volume=quantity(60,"mL","60 mL"),stage="workup"),
    participant("etoac_extract", "EtOAc", "workup_agent", e26, volume=quantity(270,"mL","3 × 90 mL"),stage="workup"),
    participant("dry", "Na2SO4", "reagent", e26), participant("celite", "Celite", "reagent", e26),
    participant("product5", "cmp_5", "product", e26),
]
triol_steps = [
    step("prepare_sia_addition", e26, ("alkene","bh3","thf_sia"), "add BH3·DMS dropwise to 2-methyl-2-butene in THF", quantity(0,"°C","0 °C")),
    step("prepare_sia_hold", e26, ("alkene","bh3","thf_sia"), "stir", quantity(0,"°C","0 °C"), quantity(2,"h","2 h")),
    step("prepare_ipc_addition", e26, ("ipc_complex","bf3","thf_ipc"), "add BF3·Et2O dropwise", "room temperature"),
    step("prepare_ipc_hold", e26, ("ipc_complex","bf3","thf_ipc"), "stir, then filter slurry", "room temperature", quantity(2.5,"h","2.5 h")),
    step("first_hydroboration_addition", e26, ("allene4","sia","thf_allene"), "add Sia2BH solution dropwise", quantity(-78,"°C","−78 °C"), quantity(30,"min","30 min")),
    step("warm_minus50", e26, ("allene4","sia"), "warm and hold", quantity(-50,"°C","−50 °C"), quantity(1,"h","1 h")),
    step("warm_minus30", e26, ("allene4","sia"), "warm and hold", quantity(-30,"°C","−30 °C"), quantity(1,"h","1 h")),
    step("second_hydroboration", e26, ("ipc","allene4"), "add (+)-IpcBH2, warm to room temperature and stir", "room temperature", quantity(5,"h","5 h")),
    step("methanol_quench", e26, ("methanol",), "quench by MeOH addition; gas evolution caution", quantity(0,"°C","0 °C")),
    step("oxidation", e26, ("base","peroxide"), "slowly add 3 N NaOH then 30% H2O2", "room temperature", quantity(3,"h","3 h")),
    step("workup", e26, ("etoac_add","etoac_extract","celite","dry"), "add EtOAc, filter through Celite, extract, wash, dry, and chromatograph")]
reaction("scheme2_4_to_5", (2,6), "specific_experiment", ("scale_up",), triol_parts, triol_steps,
    [outcome("yield_5", "yield", ("cmp_5",), e26, quantity(72,"%","72%; 891 mg"), basis="isolated after flash chromatography", yield_type="isolated"),
     outcome("dr_5", "selectivity", ("cmp_5",), e26, ("dr >20:1", {"metric":"dr","ratio":{"major":">20","minor":1},"order":["syn,anti,syn-triol 5","other diastereomer(s)"],"method":None}))],
    transformation="stereoselective double hydroboration/oxidation of 2,3-allenol", solvent="THF",
    scale=("allene 4 (3.0 mmol)", {"number":3.0,"unit":"mmol","reference_molecule_id":"cmp_4"}),
    procedure="Full staged procedure on PDF page 6; Scheme 2 supplies the dr and product relationship.",
    unresolved=("Exact Sia2BH concentration is not explicitly stated; (+)-IpcBH2 is reported as approximately 3.5 equivalents.",
                "The dr measurement method is not stated in the reviewed pages."),
    changes=(("C","O",0,1,"formed"),("C","O",0,1,"formed"),("C","H",0,1,"hydrogen_gain")), location="Scheme 2 and page-6 triol 5 procedure")

# Scheme 3 entry 1 is the same experiment as Table 1 row 1; evidence is merged below.
e23 = ev(2,3)
reaction("table1_entry1", (2,3), "specific_experiment", ("optimization",),
    [participant("triol5","cmp_5","reactant",e23), participant("nah","NaH","reagent",e23,equivalents=quantity(4,"equiv","4 eq"),reference="cmp_5"),
     participant("tris","Trisimid","reagent",e23), participant("thf","THF","solvent",e23), participant("internal","CH2Br2","reagent",e23),
     participant("p9","cmp_9","product",e23), participant("p10","cmp_10","product",e23)],
    [step("cyclization",e23,("triol5","nah","tris","thf","p9","p10"),"sulfonylation/cyclization",quantity(-25,"°C","−25 °C"),quantity(20,"h","20 h"))],
    [outcome("conversion","conversion",("cmp_5",),e23,quantity(100,"%","100%"),basis="1H NMR (400 MHz, CDCl3)"),
     outcome("ratio_9_10","selectivity",("cmp_9","cmp_10"),e23,("9/10 = 2/1",{"metric":"product_ratio","ratio":[2,1],"order":["9","10"],"method":"1H NMR (400 MHz, CDCl3)"}),basis="1H NMR (400 MHz, CDCl3)"),
     outcome("combined_yield","yield",("cmp_9","cmp_10"),e23,quantity(77,"%","77%"),basis="1H NMR with CH2Br2 internal standard",yield_type="combined NMR yield")],
    transformation="regioselective epoxide formation", scale=("400-mg scale (footnote e)",{"raw":"400 mg","number":400,"unit":"mg","reference_molecule_id":"cmp_5"}),
    unresolved=("Individual yields of epoxides 9 and 10 are not reported for this row.",), changes=(("C","O",0,1,"formed"),), location="Scheme 3a and Table 1 entry 1")

reaction("scheme3_methylation_9_10", (2,3), "specific_experiment", ("other",),
    [participant("mix9","cmp_9","reactant",e23),participant("mix10","cmp_10","reactant",e23),
     participant("mei","MeI","reagent",e23),participant("nah","NaH","reagent",e23),participant("thf","THF","solvent",e23),
     participant("p3","cmp_3","product",e23),participant("p3a","cmp_3a","product",e23)],
    [step("methylation",e23,("mix9","mix10","mei","nah","thf","p3","p3a"),"methylation in THF")],
    [outcome("combined_yield_3_3a","yield",("cmp_3","cmp_3a"),e23,quantity(95,"%","95%"),basis=None,scope="combined outcome shown for the mixed epoxides")],
    transformation="O-methylation", unresolved=("Temperature, time, measurement basis, and individual product yields are not reported in the reviewed pages.",),
    changes=(("C","O",0,1,"formed"),), location="Scheme 3a, epoxides 9/10 to 3/3a")

scheme_record("scheme3b_sulfonyl_migration", (2,), "proposed sulfonyl migration and stereoinvertive epoxide formation",
    [participant("s5","cmp_5","reactant",e2),participant("p10","cmp_10","product",e2)],
    notes=("Authors describe this pathway as presumptive; it rationalizes formation of C7-epimer 10.",),
    changes=(("C","O",0,1,"formed"),), record_type="proposed_mechanism")


# Remaining Table 1 rows: each condition row is an independent optimization record.
table_rows = [
    # no, reagent, equiv, temperature, time, conversion, ratio, combined yield
    (2,"TsCl",4,"rt",1.5,66,"4.5/1",17),
    (3,"NsCl",4,"rt",1.5,100,"8/1",18),
    (4,"anisyl_SO2Cl",4,"rt",1.5,100,"table dash (–)","table dash (–)"),
    (5,"TrisCl",4,"rt",1.5,83,"1/2",36),
    (6,"Trisimid",1.3,"0 °C",1,100,"1/1",90),
    (7,"Trisimid",4,"−78 to −25 °C",6,100,"2/1",84),
    (8,"DEAD",None,"reflux",4,100,"4/1",87),
    (9,"DIAD",None,"130 °C",8,100,"12/1",91),
]
for no, reagent_mid, equiv, temp_raw, hours, conv, ratio, yld in table_rows:
    evidence = ev(3)
    base_mid = "KHMDS" if no == 6 else "PPh3" if no in (8,9) else "NaH"
    parts = [participant("triol5","cmp_5","reactant",evidence),
             participant("activator",reagent_mid,"reagent",evidence),
             participant("base",base_mid,"reagent",evidence,
                         equivalents=quantity(equiv,"equiv",f"{equiv} eq") if equiv is not None else None,
                         reference="cmp_5"),
             participant("solvent","toluene" if no in (8,9) else "THF","solvent",evidence),
             participant("internal","CH2Br2","reagent",evidence),
             participant("p9","cmp_9","product",evidence), participant("p10","cmp_10","product",evidence)]
    if no in (8,9):
        parts.append(participant("side","side_thf","product",evidence))
    if temp_raw == "rt": temp = "room temperature"
    elif temp_raw == "reflux": temp = "toluene reflux (exact temperature unreported)"
    elif temp_raw.startswith("−78"):
        temp = (temp_raw,{"minimum":-78,"maximum":-25,"unit":"°C","raw_text":temp_raw})
    else:
        temp = (temp_raw,{"number":float(temp_raw.split()[0]),"unit":"°C","raw_text":temp_raw})
    outs = [outcome("conversion","conversion",("cmp_5",),evidence,quantity(conv,"%",f"{conv}%"),basis="1H NMR (400 MHz, CDCl3)"),
            outcome("ratio_9_10","selectivity",("cmp_9","cmp_10"),evidence,
                    (ratio, {"metric":"product_ratio","raw":ratio,"order":["9","10"],"method":"1H NMR (400 MHz, CDCl3)"}),
                    basis="1H NMR (400 MHz, CDCl3)"),
            outcome("combined_yield","yield",("cmp_9","cmp_10"),evidence,
                    quantity(yld,"%",f"{yld}%") if isinstance(yld,(int,float)) else (yld,yld),
                    basis="1H NMR with CH2Br2 internal standard",yield_type="combined NMR yield")]
    if no in (8,9):
        outs.append(outcome("side_product_yield","yield",("side_thf",),evidence,
                            ("~6%",{"number":6,"unit":"%","qualifier":"approx","raw_text":"~6%"}),
                            basis="footnote f; method not separately stated",yield_type="approximate side-product yield"))
    if no == 9:
        outs.append(outcome("isolated_yield_9","yield",("cmp_9",),evidence,quantity(81,"%","81%"),
                            basis="isolated yield (footnote g)",yield_type="isolated"))
    unresolved_row = ["Individual yields of epoxides 9 and 10 are not reported."]
    if no == 4:
        unresolved_row.append("Table dash is retained literally; it is not interpreted as zero or no reaction.")
    if no in (8,9):
        unresolved_row.append("The five-membered side-product structure was described in ref. 16, not structurally established in the reviewed pages.")
    reaction(f"table1_entry{no}",(3,),"specific_experiment",("optimization",),parts,
             [step("cyclization",evidence,tuple(p["participant_id"] for p in parts),"epoxide-forming cyclization",temp,quantity(hours,"h",f"{hours} h"))],
             outs,transformation="optimization of epoxide formation",unresolved=tuple(unresolved_row),
             changes=(("C","O",0,1,"formed"),),location=f"Table 1 entry {no}")


# Page 3, Scheme 4: unsuccessful concept plus the implemented alternative sequence.
e3 = ev(3)
reaction("scheme4a_bisboronate_attempt",(3,),"specific_experiment",("control","other"),
    [participant("s3","cmp_3","reactant",e3),participant("bisb","lithiated_bisboronate","reagent",e3),
     participant("base","NaOH","reagent",e3),participant("peroxide","H2O2","reagent",e3)],
    [step("attempt",e3,("s3","bisb","base","peroxide"),"attempted epoxide opening/oxidative deborylation under several varied conditions")],
    [outcome("no_success","other",(),e3,("no success after several attempts","no successful formation of 12 reported"),basis="qualitative discussion")],
    transformation="attempted one-carbon elongation to beta-hydroxy aldehyde 12",
    notes=("The source aggregates several attempts that varied reagent, temperature, and solvent; individual trials are not enumerated.",),
    unresolved=("Individual conditions, run count, and recovered/decomposition products are unreported; no product is fabricated.",),
    location="Scheme 4a and adjacent discussion")

reaction("scheme4b_3_to_13",(3,),"specific_experiment",("other",),
    [participant("s3","cmp_3","reactant",e3),participant("cu","CuI","catalyst",e3),participant("vinyl","vinylMgBr","reagent",e3),
     participant("thf","THF","solvent",e3),participant("tbs","TBSOTf","reagent",e3),participant("lut","lutidine","reagent",e3),
     participant("dcm","CH2Cl2","solvent",e3),participant("p13","cmp_13","product",e3)],
    [step("vinylation",e3,("s3","cu","vinyl","thf"),"vinylmagnesium bromide/CuI epoxide opening",("−20 to 0 °C",{"minimum":-20,"maximum":0,"unit":"°C","raw_text":"−20 to 0 °C"})),
     step("silylation",e3,("tbs","lut","dcm","p13"),"TBS protection")],
    [outcome("yield_13","yield",("cmp_13",),e3,quantity(95,"%","95%"),basis=None,scope="combined one-pot vinylation/silylation")],
    transformation="vinylative epoxide opening and silylation",unresolved=("Yield measurement basis and stage time are unreported on the scheme.",),
    changes=(("C","C",0,1,"formed"),("Si","O",0,1,"formed")),location="Scheme 4b, 3 to 13")

reaction("scheme4b_13_to_14",(3,),"specific_experiment",("other",),
    [participant("s13","cmp_13","reactant",e3),participant("o3","ozone","reagent",e3),participant("pph3","PPh3","reagent",e3),
     participant("dcm","CH2Cl2","solvent",e3),participant("p14","cmp_14","product",e3)],
    [step("ozonolysis",e3,("s13","o3","dcm"),"ozonolysis"),step("reductive_workup",e3,("pph3","p14"),"PPh3 reductive workup")],
    [outcome("yield_14","yield",("cmp_14",),e3,quantity(91,"%","91%"),basis=None)],
    transformation="ozonolysis to aldehyde",unresolved=("Temperature, time, and yield basis are unreported on this scheme.",),
    changes=(("C","O",0,2,"formed"),),location="Scheme 4b, 13 to 14")

reaction("scheme4b_14_to_16",(3,),"specific_experiment",("other",),
    [participant("s14","cmp_14","reactant",e3),participant("r15","cmp_15","reagent",e3),participant("bf3","BF3_Et2O","reagent",e3),
     participant("thf","THF","solvent",e3),participant("acyl","acryloyl_chloride","reagent",e3),participant("dcm","CH2Cl2","solvent",e3),
     participant("p16","cmp_16","product",e3)],
    [step("allylation",e3,("s14","r15","bf3","thf"),"syn pentenylation",quantity(-78,"°C","−78 °C")),
     step("acrylation",e3,("acyl","dcm","p16"),"trap alcohol with acryloyl chloride","room temperature")],
    [outcome("yield_16","yield",("cmp_16",),e3,quantity(71,"%","71%"),basis=None,scope="one-pot allylation/acrylation"),
     outcome("dr_16","selectivity",("cmp_16",),e3,("13:1",{"metric":"dr","ratio":[13,1],"order":["syn stereoisomer 16","other diastereomer"],"method":None}))],
    transformation="stereoselective pentenylation and acrylation",unresolved=("Measurement basis for yield and dr is not reported in the reviewed pages.",),
    changes=(("C","C",0,1,"formed"),("C","O",0,1,"formed")),location="Scheme 4b, 14 to 16")

# 16 -> 17 procedure is printed on page 6; merge exact charges with Scheme 4.
e36 = ev(3,6)
reaction("scheme4b_16_to_17",(3,6),"specific_experiment",("scale_up",),
    [participant("s16","cmp_16","reactant",e36,mass=quantity(93,"mg","93 mg"),amount=quantity(0.136,"mmol","0.136 mmol")),
     participant("grubbs","GrubbsII","catalyst",e36,mass=quantity(23,"mg","23 mg"),amount=quantity(0.027,"mmol","0.027 mmol")),
     participant("dcm","CH2Cl2","solvent",e36,volume=quantity(28,"mL","27 mL + 1.0 mL")),participant("p17","cmp_17","product",e36)],
    [step("heat_substrate",e36,("s16","dcm"),"stir substrate solution",quantity(40,"°C","40 °C")),
     step("catalyst_addition",e36,("grubbs","dcm"),"slowly add Grubbs-II solution; monitor by TLC",quantity(40,"°C","40 °C")),
     step("purification",e36,("p17",),"remove solvent and flash chromatograph on silica, PE/EA 15:1")],
    [outcome("yield_17","yield",("cmp_17",),e36,quantity(88,"%","78 mg, 88%"),basis="isolated after flash chromatography",yield_type="isolated")],
    transformation="ring-closing metathesis to pyrone 17",solvent="freshly distilled CH2Cl2",
    scale=("ester 16, 0.136 mmol",{"number":0.136,"unit":"mmol","reference_molecule_id":"cmp_16"}),
    unresolved=("Reaction completion time is reported only as TLC monitoring, not a numeric duration.",),
    changes=(("C","C",0,2,"formed"),),location="Scheme 4b and page-6 compound 17 procedure")


def simple_arrow(suffix, pages, source_mid, product_mid, reagents, yld=None, temp=None, time=None,
                 transformation=None, notes=(), unresolved=(), changes=(), categories=("other",), basis=None,
                 location=None):
    evidence = ev(*pages)
    parts = [participant("substrate",source_mid,"reactant",evidence)]
    solvent_mids = {"THF","CH2Cl2","MeOH","toluene","EtOH","hexane","Et2O","EtOAc"}
    for i, mid in enumerate(reagents,1):
        parts.append(participant(f"agent{i}",mid,"solvent" if mid in solvent_mids else "reagent",evidence))
    if product_mid:
        parts.append(participant("product",product_mid,"product",evidence))
    outs = []
    if yld is not None and product_mid:
        outs.append(outcome("yield","yield",(product_mid,),evidence,
                            quantity(yld,"%",f"{yld}%") if isinstance(yld,(int,float)) else yld,
                            basis=basis))
    return reaction(suffix,pages,"specific_experiment",categories,parts,
                    [step("reaction",evidence,tuple(p["participant_id"] for p in parts),transformation,temp,time)],
                    outs,transformation=transformation,notes=notes,
                    unresolved=tuple(unresolved) + (() if basis is not None or yld is None else ("Yield measurement basis is unreported in the reviewed pages.",)),
                    changes=changes,location=location)


# Page 4, Schemes 5-7.
simple_arrow("scheme5_17_to_18",(3,4),"cmp_17","cmp_18",("TBAF","THF"),75,
             transformation="TBAF-promoted pyrone fragmentation",
             notes=("Desired mono-desilylated alcohol 19 was not obtained.",),
             unresolved=("The 75% isolated procedure appears later in the article; page 4 supplies the reaction identity.",),
             changes=(("C","O",1,0,"broken"),),location="Scheme 5, 17 to 18")

simple_arrow("scheme5_17_to_20",(3,4),"cmp_17","cmp_20",("TsOH","MeOH"),80,
             temp=quantity(60,"°C","60 °C"),transformation="acidic removal of two silyl groups",
             notes=("No desired compound 19 was detected.",),
             unresolved=("Reaction time and yield measurement basis are not stated in the reviewed pages.",),
             changes=(),location="Scheme 5, 17 to 20")

simple_arrow("scheme5_20_to_21",(3,4),"cmp_20","cmp_21",("TEMPO","BAIB","CH2Cl2"),79,
             temp="room temperature",transformation="oxidation/lactonization to bislactone",
             notes=("Lactol 20a is proposed as an intermediate, not isolated as a standalone experiment.",),
             unresolved=("Reaction time and yield measurement basis are not reported in pages 1-6.",),changes=(("C","O",0,2,"formed"),("C","O",0,1,"formed")),
             location="Scheme 5, 20 via 20a to 21")

simple_arrow("scheme6_3a_to_22",(4,),"cmp_3a","cmp_22",("CuI","vinylMgBr","THF","lutidine","TBSOTf"),95,
             temp=("−20 to 0 °C",{"minimum":-20,"maximum":0,"unit":"°C","raw_text":"−20 to 0 °C"}),
             transformation="vinylative epoxide opening and silylation of C7 epimer",changes=(("C","C",0,1,"formed"),("Si","O",0,1,"formed")),
             location="Scheme 6, 3a to 22")
simple_arrow("scheme6_22_to_23",(3,4),"cmp_22","cmp_23",("ozone","PPh3","CH2Cl2"),87,
             transformation="ozonolysis to aldehyde 23",changes=(("C","O",0,2,"formed"),),location="Scheme 6, 22 to 23")
simple_arrow("scheme6_23_to_24",(3,4),"cmp_23","cmp_24",("cmp_15","BF3_Et2O","THF","acryloyl_chloride","CH2Cl2"),65,
             temp=quantity(-78,"°C","−78 °C"),transformation="pentenylation/acrylation of C7 epimer",
             changes=(("C","C",0,1,"formed"),("C","O",0,1,"formed")),location="Scheme 6, 23 to 24")
simple_arrow("scheme6_24_to_25",(4,),"cmp_24","cmp_25",("GrubbsII","CH2Cl2"),80,
             temp=quantity(40,"°C","40 °C"),transformation="ring-closing metathesis to pyrone 25",
             notes=("91% brsm is also reported.",),changes=(("C","C",0,2,"formed"),),location="Scheme 6, 24 to 25")

e4 = ev(4)
reaction("scheme6_25_desilylation_branches",(4,),"specific_experiment",("other",),
    [participant("s25","cmp_25","reactant",e4),participant("tbaf","TBAF","reagent",e4),participant("thf","THF","solvent",e4),
     participant("p26a","cmp_26a","product",e4),participant("p26b","cmp_26b","product",e4),participant("p26c","cmp_26c","product",e4)],
    [step("desilylation",e4,("s25","tbaf","thf","p26a","p26b","p26c"),"TBAF treatment in THF")],
    [outcome("yield_26a","yield",("cmp_26a",),e4,quantity(37,"%","37%"),basis=None),
     outcome("yield_26b","yield",("cmp_26b",),e4,quantity(20,"%","20%"),basis=None),
     outcome("yield_26c","yield",("cmp_26c",),e4,quantity(32,"%","32%"),basis=None)],
    transformation="TBAF desilylation with fragmentation and oxa-Michael branches",
    notes=("26c is attributed to oxa-Michael cyclization; this mechanistic attribution is not a separate experiment.",),
    unresolved=("Exact charges, time, temperature, yield measurement basis, distribution accounting, and unisolated remainder are not reported on page 4.",),
    changes=(("C","O",0,1,"formed"),),location="Scheme 6, 25 to 26a/26b/26c")

# Scheme 7 distinguishes one observed tosylation from three failed desired transformations.
simple_arrow("scheme7_26a_ph3pi2_failed",(4,),"cmp_26a",None,("Ph3PI2",),None,
             transformation="attempted iodination to 27",notes=("Desired iodide 27 was not produced.",),
             unresolved=("Outcome products, time, temperature, and exact loadings are unreported.",),categories=("control","other"),
             location="Scheme 7, Ph3P-I2 attempt")
simple_arrow("scheme7_26a_to_27a",(4,),"cmp_26a","cmp_27a",("TsCl","Et3N"),None,
             transformation="tosylation",unresolved=("Yield, time, temperature, and measurement basis are unreported.",),
             changes=(("S","O",0,1,"formed"),),location="Scheme 7, 26a to 27a")
simple_arrow("scheme7_27a_nai_failed",(4,),"cmp_27a",None,("NaI","acetone"),None,temp="acetone reflux (exact temperature unreported)",
             transformation="attempted iodide substitution to 27",notes=("Desired iodide 27 was not produced.",),
             unresolved=("Exact loadings/time and outcome products are unreported.",),categories=("control","other"),
             location="Scheme 7, NaI attempt")
simple_arrow("scheme7_27a_crosscoupling_failed",(4,),"cmp_27a",None,(),None,
             transformation="attempted transition-metal-catalyzed cross-coupling",
             notes=("Various transition-metal-catalyzed conditions failed to install the C12-C14 vinylic motif.",),
             unresolved=("Catalysts, partners, individual trial conditions, and products are not enumerated in the reviewed pages.",),categories=("control","other"),
             location="Scheme 7 and adjacent discussion")


# Page 5, Schemes 8-10. Same-product arrows under materially different conditions remain separate.
e5 = ev(5)
simple_arrow("scheme8_16_to_28",(4,5),"cmp_16","cmp_28",("TBAF","THF"),86,
             transformation="selective TBDPS deprotection",notes=("86% is reported brsm.",),
             unresolved=("Temperature, time, and unrecovered isolated yield are not reported in the reviewed pages.",),
             location="Scheme 8, 16 to 28")

reaction("scheme8_28_to_29",(4,5),"specific_experiment",("other",),
    [participant("s28","cmp_28","reactant",ev(4,5)),participant("dmp","DMP","reagent",ev(4,5)),
     participant("dcm","CH2Cl2","solvent",ev(4,5)),participant("phos","phosphonate","reagent",ev(4,5)),
     participant("base","tBuOK","reagent",ev(4,5)),participant("thf","THF","solvent",ev(4,5)),participant("p29","cmp_29","product",ev(4,5))],
    [step("oxidation",ev(4,5),("s28","dmp","dcm"),"Dess-Martin oxidation"),
     step("hwe",ev(4,5),("phos","base","thf","p29"),"Horner-Wadsworth-Emmons olefination")],
    [outcome("direct_yield","yield",("cmp_29",),ev(4,5),quantity(47,"%","47%"),basis=None,scope="single pass"),
     outcome("recycle_yield","yield",("cmp_29",),ev(4,5),quantity(80,"%","80% after two cycles"),basis=None,scope="combined after two recycle cycles")],
    transformation="oxidation/HWE homologation to E-enone 29",
    notes=("Phosphate Michael addition to the acrylate consumed desired 29; recovered starting material was recycled.",),
    unresolved=("Scheme also shows 86% brsm for the single pass; isolated-versus-brsm accounting and measurement bases require review.",
                "The Michael-adduct structure is discussed but not registered as a source-identified product in pages 1-6."),
    changes=(("C","C",0,2,"formed"),),location="Scheme 8, 28 to 29 and adjacent discussion")

simple_arrow("scheme8_29_to_30",(5,),"cmp_29","cmp_30",("GrubbsII","CH2Cl2"),53,
             transformation="ring-closing metathesis",notes=("93% brsm is also shown.",),
             unresolved=("Temperature, time, and yield basis are not stated on page 5.",),
             changes=(("C","C",0,2,"formed"),),location="Scheme 8, 29 to 30")

simple_arrow("scheme8_30_to_31",(4,5),"cmp_30","cmp_31",("TsNHNH2","MeOH","NaBH4","AcOH"),85,
             temp=quantity(60,"°C","60 °C for hydrazone formation"),transformation="hydrazone reduction/diazene rearrangement cascade",
             unresolved=("Detailed stage times and yield measurement basis are not reported in pages 1-6.",),
             changes=(("C","H",0,1,"hydrogen_gain"),),location="Scheme 8, 30 to 31")

reaction("scheme8_30_to_pironetin_onepot",(4,5),"specific_experiment",("other",),
    [participant("s30","cmp_30","reactant",ev(4,5)),participant("hydrazide","TsNHNH2","reagent",ev(4,5)),
     participant("meoh","MeOH","solvent",ev(4,5)),participant("borohydride","NaBH4","reagent",ev(4,5)),
     participant("acoh","AcOH","solvent",ev(4,5)),participant("hcl","HCl","reagent",ev(4,5)),
     participant("ethanol","EtOH","solvent",ev(4,5)),participant("intermediate","cmp_31","unknown",ev(4,5)),
     participant("p1","cmp_1","product",ev(4,5))],
    [step("reduction_rearrangement",ev(4,5),("s30","hydrazide","meoh","borohydride","acoh","intermediate"),"hydrazone formation, reduction, and diazene rearrangement"),
     step("desilylation",ev(4,5),("intermediate","hcl","ethanol","p1"),"acidic TBS removal")],
    [outcome("overall_yield","yield",("cmp_1",),ev(4,5),quantity(76,"%","76% over two operations"),basis=None,
             scope="overall yield from enone-pyrone 30 to pironetin 1")],
    transformation="one-pot enone reduction/rearrangement and desilylation",
    unresolved=("The 76% is retained as a start-to-end overall yield and is not assigned to the final deprotection alone.",
                "Exact operation-level yield basis is unreported in pages 1-6."),
    changes=(("C","H",0,1,"hydrogen_gain"),),location="Scheme 8, 30 to pironetin 1")

simple_arrow("scheme9_26a_to_32",(4,5),"cmp_26a","cmp_32",("DMP","CH2Cl2","phosphonate","tBuOK","THF"),84,
             transformation="oxidation/HWE homologation of C7 epimer",changes=(("C","C",0,2,"formed"),),
             location="Scheme 9, 26a to 32")
simple_arrow("scheme9_32_to_33",(5,),"cmp_32","cmp_33",("TsNHNH2","MeOH","NaBH4","AcOH"),None,
             transformation="hydrazone reduction/diazene rearrangement cascade",
             unresolved=("The page states that the complete carbon framework was established, but no isolated yield is given in the reviewed pages.",),
             changes=(("C","H",0,1,"hydrogen_gain"),),location="Scheme 9, 32 to 33")
simple_arrow("scheme9_33_acid_failure",(5,),"cmp_33",None,("HCl","EtOH"),None,
             transformation="attempted acidic TBS removal to 7-epi-pironetin",
             notes=("Several unidentified products formed; desired 7-epi-pironetin was not reported.",),
             unresolved=("Unidentified products, their amounts, and conversion are not reported.",),categories=("control","other"),
             location="Scheme 9, failed deprotection of 33")

simple_arrow("scheme10_14_crotonoyl_route",(5,),"cmp_14","cmp_34",("cmp_15","BF3_Et2O","THF","crotonoyl_chloride","Et3N","CH2Cl2"),None,
             temp=quantity(-78,"°C","−78 °C for allylation; then room temperature"),
             transformation="pentenylation/crotonoylation with double-bond migration",
             notes=("34 was inseparable from the initially crotonylated derivative; screening Et3N and DIPEA did not efficiently prevent isomerization.",),
             unresolved=("No isolated yield, mixture ratio, time, or measurement method is reported on page 5.",),
             changes=(("C","C",0,1,"formed"),("C","O",0,1,"formed")),location="Scheme 10, 14 to 34")

reaction("scheme10_14_to_35",(5,),"specific_experiment",("other",),
    [participant("s14","cmp_14","reactant",e5),participant("r15","cmp_15","reagent",e5),participant("bf3","BF3_Et2O","reagent",e5),
     participant("thf","THF","solvent",e5),participant("cinn","cinnamoyl_chloride","reagent",e5),participant("et3n","Et3N","reagent",e5),
     participant("dcm","CH2Cl2","solvent",e5),participant("p35","cmp_35","product",e5)],
    [step("allylation",e5,("s14","r15","bf3","thf"),"pentenylation",quantity(-78,"°C","−78 °C")),
     step("cinnamoylation",e5,("cinn","et3n","dcm","p35"),"cinnamoylation","room temperature")],
    [outcome("yield_35","yield",("cmp_35",),e5,basis=None,
             conflict_alternatives=(("68% (Scheme 10)",{"number":68,"unit":"%","raw_text":"68%"},e5),
                                    ("78% (adjacent prose)",{"number":78,"unit":"%","raw_text":"78%"},e5)))],
    transformation="stereoselective pentenylation/cinnamoylation",
    notes=("The scheme reports 68%, while the adjacent prose reports 78%; both page-5 values are retained without choosing.",),
    unresolved=("Source conflict: 68% versus 78% overall yield; measurement basis is not reported in the reviewed pages.",),
    changes=(("C","C",0,1,"formed"),("C","O",0,1,"formed")),location="Scheme 10 and adjacent prose, 14 to 35")

simple_arrow("scheme10_35_to_36",(5,),"cmp_35","cmp_36",("TBAF","THF"),82,
             transformation="selective TBDPS deprotection",
             notes=("No appreciable diol contamination was reported.",),location="Scheme 10, 35 to 36")
simple_arrow("scheme10_36_to_37",(5,),"cmp_36","cmp_37",("DMP","CH2Cl2","phosphonate","tBuOK","THF"),78,
             temp=quantity(45,"°C","45 °C for the HWE stage"),transformation="oxidation/HWE homologation to E-enone 37",
             notes=("No denoted phosphate Michael-addition impurity was observed.",),basis="isolated yield",
             changes=(("C","C",0,2,"formed"),),location="Scheme 10, 36 to 37")
simple_arrow("scheme10_37_to_30",(5,),"cmp_37","cmp_30",("GrubbsII","toluene"),80,
             temp=quantity(100,"°C","100 °C"),transformation="ring-closing metathesis to enone-pyrone 30",
             basis="isolated yield",changes=(("C","C",0,2,"formed"),),location="Scheme 10, 37 to 30")


# Page 6 Mosher ester preparations are two experiments because reagent configuration and isolated product differ.
def mosher_record(configuration, acid_mid, product_mid, mass, pct):
    evidence = ev(2,6)
    parts = [
        participant("alcohol","alcohol_9","reactant",evidence,mass=quantity(6.4,"mg","6.4 mg"),amount=quantity(0.014,"mmol","0.014 mmol")),
        participant("acid",acid_mid,"reagent",evidence,mass=quantity(17,"mg","17 mg"),amount=quantity(0.073,"mmol","0.073 mmol")),
        participant("dmf","DMF","reagent",evidence,volume=quantity(2.5,"uL","2.5 μL"),amount=quantity(0.03,"mmol","0.03 mmol")),
        participant("oxalyl","oxalyl_chloride","reagent",evidence,volume=quantity(31,"uL","31 μL"),amount=quantity(0.36,"mmol","0.36 mmol")),
        participant("hexane","hexane","solvent",evidence,volume=quantity(0.7,"mL","0.7 mL")),
        participant("dcm","CH2Cl2","solvent",evidence,volume=quantity(2.0,"mL","1.0 mL + 1.0 mL")),
        participant("et3n","Et3N","reagent",evidence,volume=quantity(21,"uL","21 μL"),amount=quantity(0.14,"mmol","0.14 mmol")),
        participant("dmap","DMAP","catalyst",evidence,mass=quantity(1.0,"mg","1.0 mg")),
        participant("ether","Et2O","workup_agent",evidence,volume=quantity(5,"mL","5 mL"),stage="workup"),
        participant("water","water","workup_agent",evidence,stage="workup"),participant("brine","brine","workup_agent",evidence,stage="workup"),
        participant("dry","Na2SO4","workup_agent",evidence,stage="workup"),participant("pe","petroleum_ether","purification_agent",evidence,stage="purification"),
        participant("ethylacetate","EtOAc","purification_agent",evidence,stage="purification"),participant("product",product_mid,"product",evidence)]
    return reaction(f"page6_{configuration.lower()}_mtpa_11",(2,6),"specific_experiment",("derivatization",),parts,
        [step("acid_chloride_formation",evidence,("acid","dmf","oxalyl","hexane"),"stir MTPA acid with catalytic DMF and oxalyl chloride",None,quantity(2,"h","2 h")),
         step("filter_concentrate",evidence,("acid","hexane"),"filter and concentrate to crude acid chloride"),
         step("esterification",evidence,("alcohol","acid","et3n","dmap","dcm","product"),"add acid chloride solution to alcohol/Et3N/DMAP solution and stir",None,quantity(6,"h","6 h")),
         step("workup",evidence,("ether","water","brine","dry"),"add ether; wash with water and brine; dry over Na2SO4",stage="workup"),
         step("purification",evidence,("pe","ethylacetate","product"),"PTLC, PE/EA 30:1",stage="purification")],
        [outcome("isolated_yield","yield",(product_mid,),evidence,quantity(pct,"%",f"{mass} mg, {pct}%"),
                 basis="isolated after PTLC",yield_type="isolated")],
        transformation=f"formation of the {configuration}-MTPA ester for Mosher analysis",
        scale=("alcohol derived from major epoxide 9, 0.014 mmol",{"number":0.014,"unit":"mmol","reference_molecule_id":"alcohol_9"}),
        procedure="Two-stage acid chloride preparation and esterification printed on PDF page 6.",
        unresolved=("The preparation of the alcohol from major epoxide 9 is referenced but not described on pages 1-6.",
                    "Mosher interpretation is stated to be in online Supporting Information, which was not supplied and is not associated."),
        changes=(("C","O",0,1,"formed"),),location=f"Page 6, {configuration}-MTPA-11 procedure")


mosher_record("R","R_MTPA","cmp_R11",5.2,55)
mosher_record("S","S_MTPA","cmp_S11",6.0,63)


# Non-synthetic and characterization material is retained separately from reaction experiments.
bundle["source_context"] = [
    {"context_id":"page1_biological_context","document_id":DOC,"pages":[1],"evidence_ids":ev(1),
     "kind":"background","summary":"Title, abstract, pironetin biological context, and prior total-synthesis discussion; no current-paper synthetic arrow or experiment table appears on page 1."},
    {"context_id":"page2_figure1_biology","document_id":DOC,"pages":[2],"evidence_ids":ev(2),
     "kind":"biological_figure","summary":"Figure 1 depicts pironetin binding to alpha-tubulin/Cys-316; it is not extracted as a synthetic experiment."},
    {"context_id":"page6_analytical_general","document_id":DOC,"pages":[6],"evidence_ids":ev(6),
     "kind":"general_and_characterization","summary":"General atmosphere, chromatography, NMR, optical rotation, HRMS, and FT-IR methods plus spectral characterizations are source context, not standalone synthesis records."},
    {"context_id":"page6_prior_work_and_si_reference","document_id":DOC,"pages":[6],"evidence_ids":ev(6),
     "kind":"source_boundary","summary":"Procedures for compounds 3, 4, 6, 7, 9, 13, 14, 16, 28, and 29 are assigned to prior work ref. 16; Mosher structural details are said to be online SI. No SI file was supplied or associated."},
]


page_regions = {
    1:[{"region":"title, abstract, introduction and page-wide body text","finding":"background/biological and prior-synthesis context; no current-paper synthetic experiment arrow or table","evidence_ids":ev(1)}],
    2:[{"region":"Figure 1","finding":"biological mode-of-action figure, not a synthetic experiment","evidence_ids":ev(2)},
       {"region":"Scheme 1","finding":"generic initial synthetic plan","evidence_ids":ev(2)},
       {"region":"Schemes 2 and 3","finding":"early route experiments, optimization-entry duplicate, Mosher derivatization link, and proposed sulfonyl-migration mechanism","evidence_ids":ev(2)}],
    3:[{"region":"Table 1 including footnotes a-g","finding":"nine independent epoxide-formation optimization rows with conversion, product ratio, NMR yield, side product, and isolated yield fields","evidence_ids":ev(3)},
       {"region":"Scheme 4 and adjacent text","finding":"aggregated unsuccessful bisboronate approach and implemented vinylation/ozonolysis/allylation/RCM sequence","evidence_ids":ev(3)}],
    4:[{"region":"Schemes 5 and 6","finding":"pyrone deprotection/fragmentation, oxidation, C7-epimer route, and three-product TBAF branch","evidence_ids":ev(4)},
       {"region":"Scheme 7 and adjacent text","finding":"observed tosylation and failed iodination/cross-coupling attempts kept separate","evidence_ids":ev(4)}],
    5:[{"region":"Schemes 8 and 9","finding":"route reordering, recycle-dependent HWE outcome, RCM, rearrangement, overall pironetin yield, and failed 7-epi deprotection","evidence_ids":ev(5)},
       {"region":"Scheme 10 and adjacent text","finding":"crotonyl isomerization route and improved cinnamate sequence; 35 yield conflict retained as 68% versus 78%","evidence_ids":ev(5)}],
    6:[{"region":"general experimental section","finding":"general atmosphere/analytical context; not converted into independent synthesis records","evidence_ids":ev(6)},
       {"region":"R- and S-MTPA-11 procedures","finding":"two separate derivatization experiments with all printed charges, staged operations, workup, and PTLC","evidence_ids":ev(6)},
       {"region":"triol 5 and compound 17 procedures","finding":"full staged double-hydroboration procedure and full RCM charges merged with their scheme records","evidence_ids":ev(6)},
       {"region":"spectral characterization blocks","finding":"characterization/source context, not standalone synthesis experiments","evidence_ids":ev(6)}],
}

page_unresolved = {
    1:["No reaction record is assigned because the page contains background rather than a current-paper experimental arrow."],
    2:["Several scheme-only yields lack explicit measurement bases; dense structures remain label-linked without unsafe full SMILES."],
    3:["Table dash is not interpreted; individual unsuccessful bisboronate runs are not enumerated by the source."],
    4:["Failed Scheme 7 trials omit individual catalysts/loadings/products; these remain unresolved rather than guessed."],
    5:["Compound 35 yield conflicts between Scheme 10 (68%) and prose (78%); measurement bases are often unreported."],
    6:["Online SI was not supplied; prior-work procedures and detailed Mosher structural determination were not imported."],
}

coverage_pages = []
for page in range(1,7):
    ids = [r["reaction_id"] for r in bundle["reactions"] if page in r["location"]["normalized_value"].get("pages", [])]
    coverage_pages.append({"page":page,"state":"attempted","regions":page_regions[page],"record_ids":ids,
                           "unresolved":page_unresolved[page]})
coverage = {"document_id":DOC,"semantic_accuracy_verified":False,"pages":coverage_pages}

inventory = {"document_id":DOC,"pages_reviewed":[1,2,3,4,5,6],
             "reaction_record_ids":[r["reaction_id"] for r in bundle["reactions"]],
             "molecule_ids":[m["molecule_id"] for m in bundle["molecules"]],
             "source_context_ids":[x["context_id"] for x in bundle["source_context"]],
             "limitations":["Machine extraction remains needs_review.","No Supporting Information was supplied.",
                            "Dense source drawings generally remain unresolved representations with descriptive stereochemistry."]}

st.write(RUN / "submitted.json", bundle)
st.write(RUN / "raw_model/stereo_checks.json", stereo_checks)
st.write(RUN / "raw_model/page_inventory.json", inventory)
st.write(RUN / "raw_model/source_mentions.json", st.source_mentions(DATASET,["compound 35","Table 1","Supporting Information"],context=180))
st.write(RUN / "coverage/pages_0001_0006.json", coverage)
print(json.dumps({"molecules":len(bundle["molecules"]),"reactions":len(bundle["reactions"]),
                  "coverage_pages":len(coverage_pages)},ensure_ascii=False))
