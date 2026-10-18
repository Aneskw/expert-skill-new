#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for Chem. Eur. J. 2020, 26, 8983.

The script deliberately leaves structures and experimental fields unresolved when
the supplied main-text pages do not establish them. It does not read prior records.
"""
from pathlib import Path
import copy
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import standard_records as st

DOC = "doc_865f887b40814098"
RUN = ROOT / "runs/current"
TEMPLATE = RUN / "raw_model/standard_init_page_0007.json"
bundle = json.loads(TEMPLATE.read_text())


def ev(page, image=True):
    return f"{DOC}:page{page}:{'image' if image else 'text'}"


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation, status=status)


def inferred(value, evidence, derivation):
    return f(None, value, evidence, method="inference", derivation=derivation)


def set_name(m, label, name, evidence):
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(label, label, evidence)
    if name:
        m["names"]["english"] = f(name, name, evidence)
    m["identity_note"] = inferred(
        "Machine transcription from the supplied main-text drawing or label; identity and structure require expert review.",
        evidence,
        "The source label/drawing was inspected, but no Supporting Information was supplied.",
    )


molecules = {}


def mol(key, label=None, name=None, evidence=(ev(1),), smiles=None,
        chirality="unresolved", geometry="unresolved", representation=None,
        stereo_detail=None):
    if key in molecules:
        return molecules[key]["molecule_id"]
    mid = f"{DOC}:mol:{key}"
    m = st.new_molecule(mid, DOC, evidence)
    set_name(m, label or key, name, evidence)
    if smiles:
        m["identity_status"] = "source_identified"
        m["representation_type"] = "smiles"
        m["smiles"] = f(
            {"source_structure": f"drawing/label {label or key}"}, smiles, evidence,
            method="structure_derivation",
            derivation="Candidate SMILES manually transcribed from the supplied source drawing/standard labeled small molecule; parsing is not source verification.",
        )
        representation = representation or "partial"
        if representation == "encoded":
            m["stereochemistry"]["isomeric_smiles"] = f(
                {"source_structure": f"drawing/label {label or key}"}, smiles, evidence,
                method="structure_derivation",
                derivation="The candidate SMILES encodes the source-shown alkene geometry.",
            )
    m["stereochemistry"]["chirality_status"] = inferred(
        chirality, evidence,
        "Status reflects the source drawing and the limits of the candidate representation; absence of a label was not treated as racemic.",
    )
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, evidence,
        "Status reflects whether E/Z or another geometry is explicitly shown/reported for this source identity.",
    )
    m["stereochemistry"]["representation_status"] = inferred(
        representation or ("partial" if smiles else "unresolved"), evidence,
        "Representation coverage was assessed against the supplied main-text drawing; unresolved substituents/stereochemistry were not invented.",
    )
    if stereo_detail:
        for field, value in stereo_detail.items():
            m["stereochemistry"][field] = f(value, value, evidence)
    molecules[key] = m
    return mid


# Shared compounds and small molecules. Only defensible concrete structures are encoded.
M1 = mol("1", "1", "2-(2-fluoropyridinium-1-yl)-1,1-bis(triflyl)ethan-1-ide (Yanai's reagent)",
         (ev(1), ev(2), ev(3), ev(4)))
MECHN = mol("MeCN", "CH3CN", "acetonitrile", (ev(2), ev(3), ev(4)), "CC#N",
            chirality="achiral", geometry="not_applicable", representation="not_applicable")
M4A = mol("4a", "4a", "allenol 4a", (ev(2), ev(4), ev(5)),
          "OC(c1ccccc1)C=C=C", chirality="unspecified", geometry="not_applicable",
          representation="partial",
          stereo_detail={"atom_centers": [{"site": "carbinol carbon", "configuration": "not assigned in the supplied pages"}]})
M5A = mol("5a", "5a", "sodium salt of (E)-bis(triflyl)enone 5a", (ev(2), ev(3), ev(4), ev(5)),
          "c1ccccc1/C=C/C(=O)C[C-](S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F.[Na+]",
          chirality="achiral", geometry="specified", representation="encoded",
          stereo_detail={"double_bonds": [{"alkene": "enone C=C", "configuration": "E"}],
                         "assignment": "Source text calls 5a the (E)-isomer and states that the (Z)-isomer was not obtained."})
MTFALKENE = mol("Tf2CCH2", "Tf2C=CH2", "1,1-bis(triflyl)ethene", (ev(1), ev(5)))
MWATER = mol("water", "H2O", "water", (ev(4), ev(5)), "O", chirality="achiral",
             geometry="not_applicable", representation="not_applicable")


def ensure_label(key, page, name=None, **kw):
    return mol(key, key, name or key, (ev(page),), **kw)


def participant(pid, mid, role, evidence, stage="reaction", equiv=None, molpct=None, volume=None):
    p = st.new_participant(pid, mid, role, evidence, stage)
    if equiv is not None:
        p["charges"]["equivalents"] = f(str(equiv) + " equiv", {"number": equiv, "unit": "equiv"}, evidence)
        p["equivalence_reference"] = f("relative to limiting substrate", next(
            (x for x in [M4A] if x), M4A), evidence)
    if molpct is not None:
        p["charges"]["mol_percent"] = f(str(molpct) + " mol%", {"number": molpct, "unit": "mol%"}, evidence)
        p["equivalence_reference"] = f("relative to limiting substrate", M4A, evidence)
    if volume is not None:
        p["charges"]["volume"] = f(str(volume[0]) + " " + volume[1], {"number": volume[0], "unit": volume[1]}, evidence)
    return p


def step(step_id, evidence, participants, temperature=None, time=None, operation=None, stage="reaction",
         atmosphere=None):
    s = st.new_step(step_id, evidence, stage)
    s["participant_ids"] = list(participants)
    if operation is not None:
        s["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        norm = temperature if isinstance(temperature, dict) else temperature
        s["temperature"] = f(str(temperature), norm, evidence)
    if time is not None:
        norm = time if isinstance(time, dict) else time
        s["time"] = f(str(time), norm, evidence)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, evidence)
    return s


def out(outcome_id, kind, targets, evidence, raw, normalized, product_role=None,
        yield_type=None, measurement=None, scope=None):
    o = st.new_outcome(outcome_id, kind, targets, evidence)
    o["value"] = f(raw, normalized, evidence)
    if product_role:
        o["product_role"] = f(product_role, product_role, evidence)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    if measurement:
        o["measurement_basis"] = f(measurement, measurement, evidence)
    if scope:
        o["scope"] = f(scope, scope, evidence)
    return o


def classification(labels, evidence, coverage="partial", limitation=None):
    c = copy.deepcopy(st.new_reaction("x", DOC, evidence)["classification"])
    c["mapped_reaction_smiles"] = f(evidence=evidence)
    if coverage == "not_applicable":
        c["bond_changes"] = []
        c["labels"] = inferred([], evidence, "No observed chemical transformation; no formation label assigned.")
        c["coverage"] = inferred("not_applicable", evidence, "The source reports no reaction/nonproductive outcome.")
    elif coverage == "unresolved":
        c["bond_changes"] = []
        c["labels"] = f(evidence=evidence)
        c["coverage"] = inferred("unresolved", evidence, "The arrow is a project aim or lacks a defensible observed atom mapping.")
    else:
        changes = []
        label_map = {}
        for idx, code in enumerate(labels, 1):
            prefix = "DB-" if code.startswith("DB-") else "TB-" if code.startswith("TB-") else ""
            pair = code.removeprefix("DB-").removeprefix("TB-").split("-")
            if code.startswith("other-"):
                pair = code.split("-")[1:]
            order = 2 if prefix == "DB-" else 3 if prefix == "TB-" else 1
            change = {"change_id": f"b{idx}", "atoms": [{"atom_map": None, "element": pair[0]},
                      {"atom_map": None, "element": pair[1]}], "before_order": 0,
                      "after_order": order, "change_type": "formed"}
            changes.append(f(None, change, evidence, method="structure_derivation",
                             derivation="Partial net bond change inferred from the source transformation drawing; atom mapping and coproduct balance are unresolved."))
            lab = st.bond_label(pair[0], pair[1], order)
            label_map.setdefault(lab["code"], dict(lab, supporting_change_ids=[]))["supporting_change_ids"].append(f"b{idx}")
        c["bond_changes"] = changes
        c["labels"] = f(None, list(label_map.values()), evidence, method="structure_derivation",
                        derivation="Taxonomy labels were inferred from the listed partial bond changes, not reported by the authors.")
        c["coverage"] = inferred(coverage, evidence, "Only source-visible net bond formations are represented; no complete atom mapping was asserted.")
    c["limitations"] = [limitation or "Partial source-drawing classification; full atom mapping, broken bonds, hydrogen sources and coproducts are unresolved."]
    c["model_votes"] = []
    return c


reactions = []


def reaction(suffix, page, record_type, categories, title, transformation, parts, steps=(), outcomes=(),
             labels=(), class_coverage="partial", evidence=None, unresolved=None, general=()):
    evidence = tuple(evidence or (ev(page),))
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence)
    r["location"] = f(f"PDF page {page}, {title}", {"page": page, "region": title}, evidence)
    r["general_scheme_ids"] = list(general)
    r["reaction_name"] = f(title, title, evidence)
    r["transformation_type"] = inferred(transformation, evidence, "Transformation type inferred from the source arrow and accompanying prose.")
    r["participants"] = list(parts)
    r["steps"] = list(steps)
    solvents = [p["molecule_id"] for p in parts if p["role"] == "solvent"]
    if solvents:
        r["solvent_conditions"]["composition"] = f("source-shown solvent", solvents, evidence)
    r["outcomes"] = list(outcomes)
    r["classification"] = classification(labels, evidence, class_coverage)
    r["procedure"]["original_text"] = f(
        f"Only the conditions and outcomes printed in {title} and its linked prose are available in the supplied main text.",
        f"Only the conditions and outcomes printed in {title} and its linked prose are available in the supplied main text.", evidence)
    r["notes"] = []
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": list(unresolved or ["No explicitly associated Supporting Information was supplied; charges, workup, purification and characterization may be missing."])}
    r["completeness"] = "partial"
    reactions.append(r)
    return r


def pset(substrate_mid, product_mid, evidence, reagent_mid=M1, solvent_mid=MECHN):
    return [participant("substrate", substrate_mid, "reactant", evidence),
            participant("reagent_1", reagent_mid, "reagent", evidence),
            participant("product", product_mid, "product", evidence),
            participant("solvent", solvent_mid, "solvent", evidence)]


def cascade(suffix, page, substrate_key, product_key, yield_pct=None, selectivity=None, time_value=None,
            temperature="rt", footnote=None, outcome_text=None, evidence=None):
    evidence = tuple(evidence or (ev(page),))
    sm = ensure_label(substrate_key, page)
    geometry = "specified" if selectivity is not None else "unresolved"
    stereo = {"double_bonds": [{"configuration": "E-major", "source_selectivity": selectivity}]} if selectivity is not None else None
    pm = mol(product_key, product_key, product_key, evidence, chirality="unresolved", geometry=geometry,
             representation="unresolved", stereo_detail=stereo)
    parts = pset(sm, pm, evidence)
    temp = {"number": 0, "unit": "°C"} if footnote == "0C" else ({"number": -40, "unit": "°C"} if footnote == "-40C" else "rt")
    operation = "reaction initiated at -40 °C, then allowed to proceed under the standard room-temperature conditions" if footnote == "-40C" else "reaction under the scheme conditions"
    sts = [step("reaction", evidence, [p["participant_id"] for p in parts], temp,
                {"number": time_value[0], "unit": time_value[1]} if time_value else None, operation)]
    outs = []
    if yield_pct is not None:
        outs.append(out("yield", "yield", [pm], evidence, f"{yield_pct}%", {"number": yield_pct, "unit": "%"},
                        "main product", "isolated", "isolated yield",
                        "Yield of the labeled product in the scope scheme."))
    if selectivity is not None:
        if isinstance(selectivity, (int, float)):
            raw = f"{selectivity}"
            norm = {"number": selectivity, "unit": "%"}
        elif str(selectivity).startswith(">"):
            raw = str(selectivity)
            numeric = str(selectivity)[1:].split()[0]
            norm = {"number": float(numeric), "unit": "%", "qualifier": "gt", "raw_text": str(selectivity)}
        else:
            raw = str(selectivity)
            norm = str(selectivity)
        outs.append(out("selectivity", "selectivity", [pm], evidence, raw, norm, "main product", None, None,
                        "E/Z % column; E is favored. Exact analytical method is unreported in the supplied pages."))
    if outcome_text:
        outs.append(out("reported_outcome", "other", [pm] if pm else [], evidence, outcome_text, outcome_text,
                        scope="Qualitative source outcome."))
    r = reaction(suffix, page, "specific_experiment", ["substrate_scope"], "scope entry " + product_key,
                 "trifluorosulfonylation-rearrangement cascade", parts, sts, outs,
                 ["C-C", "DB-C-O"], evidence=evidence,
                 unresolved=["Individual charges, scale, workup and purification are not given in the supplied main-text entry.",
                             "Full structures/atom mapping and all implicit-H changes remain unresolved.",
                             "The E/Z analytical method is not reported." if selectivity is not None else "No E/Z value is printed for this entry."])
    return r


# Page 1: literature context, project aim, current generic scheme and classical reactivity.
alkyne = mol("generic_alkyne", "R2-C≡C-R1", "generic alkyne", (ev(1),))
cyclobutene = mol("generic_bis_triflyl_cyclobutene", "bis(triflyl)cyclobutene", evidence=(ev(1),))
r = reaction("scheme1a_previous_literature", 1, "literature_example", [], "Scheme 1a previous literature",
             "alkyne bis(triflyl)cyclobutenation",
             [participant("alkyne", alkyne, "reactant", (ev(1),)), participant("reagent_1", M1, "reagent", (ev(1),)),
              participant("product", cyclobutene, "product", (ev(1),))], labels=["C-C", "C-C"],
             unresolved=["Cited prior work is summarized generically; conditions and yield are not reported here."])
generic_allene = mol("generic_allene", "R-substituted allene", evidence=(ev(1),))
hyp_product = mol("hypothetical_cyclobutenes", "two proposed cyclobutene regioisomers", evidence=(ev(1),),
                  chirality="unresolved", geometry="unresolved")
reaction("scheme1b_initial_aim", 1, "generic_scheme", [], "Scheme 1b initial aim of the project",
         "hypothetical allene cyclobutenation",
         [participant("allene", generic_allene, "reactant", (ev(1),)), participant("electrophile", MTFALKENE, "reagent", (ev(1),)),
          participant("hypothetical_products", hyp_product, "product", (ev(1),))],
         class_coverage="unresolved", unresolved=["This is a project aim with two hypothetical products, not a reported experiment."])
generic_allenol = mol("generic_allenol_4", "allenols 4", evidence=(ev(1), ev(5)))
generic_enone = mol("generic_enones_5", "bis(triflyl)enones 5", evidence=(ev(1), ev(5)), geometry="partially_specified",
                    stereo_detail={"assignment": "The paper reports predominantly E products for applicable scope entries, but the generic structure 5 contains variable substitution and is not a single fully assigned stereoisomer."})
reaction("scheme1c_current_generic", 1, "generic_scheme", [], "Scheme 1c current work",
         "generic trifluorosulfonylation-rearrangement cascade",
         [participant("allenol_4", generic_allenol, "reactant", (ev(1),)), participant("electrophile", MTFALKENE, "reagent", (ev(1),)),
          participant("enone_5", generic_enone, "product", (ev(1),))], labels=["C-C", "DB-C-O"],
         unresolved=["Generic R groups and atom-level changes are not fully resolved; this record carries no experimental yield."])
classical_products = mol("classical_allenol_products", "classical cyclization/rearrangement products", evidence=(ev(1),))
reaction("scheme1d_classical_reactivity", 1, "literature_example", [], "Scheme 1d classical reactivity",
         "generic electrophilic cyclization or rearrangement of allenols",
         [participant("allenol", generic_allenol, "reactant", (ev(1),)), participant("electrophile", mol("generic_Eplus", "E+", evidence=(ev(1),)), "reagent", (ev(1),)),
          participant("alternative_products", classical_products, "product", (ev(1),))], class_coverage="unresolved",
         unresolved=["The scheme summarizes cited classical reactivity and two alternative generic products; no current-paper experiment is asserted."])

# Page 2, Scheme 2: four distinct no-reaction catalyst variants, one complex mixture, one initial positive result.
M2 = mol("2", "2", "3-methylbuta-1,2-diene 2", (ev(2),), "C=C=C(C)C",
         chirality="achiral", geometry="not_applicable", representation="not_applicable")
for tag, catalyst in [("none", None), ("AgOTf", "AgOTf"), ("AuNTf2", "[(Ph3P)AuNTf2]"), ("InOTf3", "In(OTf)3")]:
    evidence = (ev(2), ev(2, False))
    parts = [participant("allene_2", M2, "reactant", evidence), participant("reagent_1", M1, "reagent", evidence),
             participant("solvent", MECHN, "solvent", evidence)]
    if catalyst:
        cm = mol("cat_" + tag, catalyst, catalyst, evidence)
        parts.append(participant("catalyst", cm, "catalyst", evidence))
    outs = [out("no_reaction", "no_reaction", [], evidence, "reaction did not proceed", True,
                scope="No reaction under this catalyst variant."),
            out("recovery", "recovery", [M2], evidence, "2 recovered", "recovered (amount unreported)",
                scope="Starting material recovery; amount and measurement basis unreported.")]
    reaction("scheme2a_no_reaction_" + tag, 2, "specific_experiment", ["control"],
             "Scheme 2a/no-reaction catalyst variant " + tag, "attempted allene bis(triflyl) functionalization",
             parts, [step("reaction", evidence, [p["participant_id"] for p in parts], "rt", None, "stirred in acetonitrile")],
             outs, class_coverage="not_applicable", evidence=evidence,
             unresolved=["Catalyst loading and reaction time are not reported; recovery quantity/method are not reported."])
M3 = ensure_label("3", 2, "allenamide 3")
parts = [participant("allenamide_3", M3, "reactant", (ev(2),)), participant("reagent_1", M1, "reagent", (ev(2),)),
         participant("solvent", MECHN, "solvent", (ev(2),))]
reaction("scheme2b_complicated_mixture", 2, "specific_experiment", ["control"], "Scheme 2b",
         "attempted allenamide bis(triflyl) functionalization", parts,
         [step("reaction", (ev(2),), [p["participant_id"] for p in parts], "rt", None, "reaction in acetonitrile")],
         [out("decomposition", "decomposition", [], (ev(2),), "complicated mixture", "complicated mixture",
              scope="No clean product identity is reported.")], class_coverage="unresolved",
         unresolved=["No clean product, yield, reaction time, scale, workup or measurement method is reported."])
parts = pset(M4A, M5A, (ev(2), ev(2, False)))
reaction("scheme2c_initial_4a_to_5a", 2, "specific_experiment", ["other"], "Scheme 2c initial experiment",
         "trifluorosulfonylation-rearrangement cascade", parts,
         [step("reaction", (ev(2),), [p["participant_id"] for p in parts], "rt", None, "reaction in acetonitrile")],
         [out("yield", "yield", [M5A], (ev(2),), "75%", {"number": 75, "unit": "%"}, "main product",
              scope="Yield of sodium salt 5a; measurement basis unreported."),
          out("selectivity", "selectivity", [M5A], (ev(2), ev(2, False)), "(Z)-isomer was not obtained",
              "Z isomer not obtained", "main product", scope="Qualitative E/Z outcome; method and detection limit unreported.")],
         labels=["C-C", "DB-C-O"], evidence=(ev(2), ev(2, False)),
         unresolved=["Reaction time, charges, scale and yield measurement basis are not reported.",
                     "The phrase 'Z-isomer was not obtained' has no analytical method or detection limit.",
                     "Product 5a was isolated as the sodium salt after silica-gel chromatography, but detailed workup is unavailable."])

# Scheme 3 scope entries and two nonproductive heterocyclic allenols.
scope3 = [
    ("4a", "5a", 80, ">99", (45, "min"), None), ("4b", "5b", 31, ">99", (5, "h"), None),
    ("4c", "5c", 45, ">99", (30, "min"), "-40C"), ("4d", "5d", 48, 97, (1, "h"), None),
    ("4e", "5e", 64, 97, (3, "h"), "-40C"), ("4f", "5f", 39, ">99", (6, "h"), None),
    ("4g", "5g", 90, None, (20, "min"), None), ("4h", "5h", 59, None, (1, "h"), None),
    ("4i", "5i", 84, None, (40, "min"), None), ("4j", "5j", 68, 88, (5, "h"), None),
]
for sk, pk, y, sel, tim, foot in scope3:
    cascade("scheme3_" + sk + "_to_" + pk, 2, sk, pk, y, sel, tim, footnote=foot)
for sk in ("5-indole", "5-oxindole"):
    sm = ensure_label(sk, 2, sk + " allenol precursor")
    evidence = (ev(2), ev(2, False))
    parts = [participant("substrate", sm, "reactant", evidence), participant("reagent_1", M1, "reagent", evidence),
             participant("solvent", MECHN, "solvent", evidence)]
    reaction("scheme3_nonproductive_" + sk.replace("-", "_"), 2, "specific_experiment", ["substrate_scope", "control"],
             "Scheme 3 nonproductive precursor " + sk, "attempted trifluorosulfonylation-rearrangement cascade", parts,
             [step("reaction", evidence, [p["participant_id"] for p in parts], "rt", None, "standard conditions applied")],
             [out("decomposition", "decomposition", [], evidence, "complicated mixtures; enone formation could not be accomplished",
                  "complicated mixture/no enone", scope="Nonproductive scope experiment.")], class_coverage="not_applicable", evidence=evidence,
             unresolved=["Time, charges, component identities of the mixture and analytical method are unreported."])

# Page 3, Scheme 4 scope and suppressed proximal-substitution cases.
scope4 = [
    ("4k", "5k", 54, 97, (1, "h"), None), ("4l", "5l", 43, ">99", (1, "h"), None),
    ("4m", "5m", 39, ">99", (7, "h"), None), ("4n", "5n", 41, ">99", (2, "h"), None),
    ("4o", "5o", 78, ">99", (2, "h"), None), ("4r", "5r", 67, ">99 (E,E)", (1, "h"), None),
    ("4s", "5s", 50, ">99", (1, "h"), None), ("4t", "5t", 58, 93, (1, "h"), None),
    ("4u", "5u", 64, 88, (1, "h"), None), ("4v", "5v", 54, ">99", (1.5, "h"), None),
    ("4w", "5w", 47, ">99", (3, "h"), "0C"), ("4x", "5x", 73, ">99", (1.5, "h"), "0C"),
]
for sk, pk, y, sel, tim, foot in scope4:
    cascade("scheme4_" + sk + "_to_" + pk, 3, sk, pk, y, sel, tim, footnote=foot)
for sk in ("4p", "4q"):
    sm = ensure_label(sk, 3, sk + " proximal-substituted allenol")
    evidence = (ev(3), ev(3, False))
    parts = [participant("substrate", sm, "reactant", evidence), participant("reagent_1", M1, "reagent", evidence),
             participant("solvent", MECHN, "solvent", evidence)]
    reaction("scheme4_suppressed_" + sk, 3, "specific_experiment", ["substrate_scope", "control"],
             "Scheme 4 enone formation suppressed for " + sk, "attempted trifluorosulfonylation-rearrangement cascade", parts,
             [step("reaction", evidence, [p["participant_id"] for p in parts], "rt", None, "standard conditions applied")],
             [out("no_reaction", "no_reaction", [], evidence, "enone formation suppressed / messy reaction", True,
                  scope="No target enone formation; mixture composition unreported.")], class_coverage="not_applicable", evidence=evidence,
             unresolved=["Time, charges, conversion, recovery, and mixture composition are not reported."])

# Scheme 5: nine separately printed deuterated product entries. The exact precursor-to-product mapping is not explicit.
deut_precursor = mol("Dn_4a_or_D2_1", "[D_n]-4a or [D2]-1", "corresponding deuterated allenol or deuterated Yanai reagent", (ev(3),),
                     chirality="unspecified", geometry="unresolved")
deut_entries = [
    ("D_5a", "[D]-5a", 61, (95, "min")), ("Dprime_5a", "[D']-5a", 62, (1, "h")),
    ("D2_5a", "[D2]-5a", 72, (1, "h")), ("D2prime_5a", "[D2']-5a", 66, (80, "min")),
    ("D3_5a", "[D3]-5a", 60, (1.5, "h")), ("D3prime_5a", "[D3']-5a", 59, (1, "h")),
    ("D3doubleprime_5a", "[D3'']-5a", 73, (70, "min")), ("D4_5a", "[D4]-5a", 63, (80, "min")),
    ("D5_5a", "[D5]-5a", 83, (75, "min")),
]
for key, label, y, tim in deut_entries:
    pm = mol(key, label, "deuterated bis(triflyl)enone " + label, (ev(3),), chirality="unresolved",
             geometry="specified", representation="unresolved",
             stereo_detail={"double_bonds": [{"configuration": "E-major", "selectivity": ">99%"}]})
    parts = [participant("deuterated_precursor", deut_precursor, "reactant", (ev(3),)),
             participant("reagent", M1, "reagent", (ev(3),)), participant("product", pm, "product", (ev(3),)),
             participant("solvent", MECHN, "solvent", (ev(3),))]
    reaction("scheme5_" + key, 3, "specific_experiment", ["substrate_scope"],
             "Scheme 5 " + label, "deuterated trifluorosulfonylation-rearrangement cascade", parts,
             [step("reaction", (ev(3),), [p["participant_id"] for p in parts], "rt", {"number": tim[0], "unit": tim[1]},
                   "reaction under standard conditions")],
             [out("yield", "yield", [pm], (ev(3),), f"{y}%", {"number": y, "unit": "%"}, "main product", "isolated", "isolated yield"),
              out("selectivity", "selectivity", [pm], (ev(3),), ">99", {"number": 99, "unit": "%", "qualifier": "gt", "raw_text": ">99"},
                  "main product", scope="E/Z % column; exact method unreported.")], labels=["C-C", "DB-C-O"],
             unresolved=["The main-text scheme groups [D_n]-4a and [D2]-1; the exact precursor identity for this individual product is not explicit.",
                         "Individual charges, deuterium-incorporation measurement details, workup and purification are unavailable."])

# Page 4, Scheme 6 transformations.
def reagent(key, label, name=None, smiles=None):
    return mol("reag_" + key, label, name or label, (ev(4),), smiles,
               chirality="achiral" if smiles else "unresolved",
               geometry="not_applicable" if smiles else "unresolved",
               representation="not_applicable" if smiles else "unresolved")

M6 = ensure_label("6", 4, "N-tosylhydrazone 6")
M7 = ensure_label("7", 4, "pyrazole 7")
M8 = ensure_label("8", 4, "allylic alcohol 8")
M9 = ensure_label("9", 4, "1,3-diene 9")
M10 = ensure_label("10", 4, "thioether 10")
M11 = ensure_label("11", 4, "pyrrole 11")
M12 = ensure_label("12", 4, "cyclobutene 12")
M13 = ensure_label("13", 4, "chromene derivative 13")
M14 = ensure_label("14", 4, "triflyl-triazole 14")
R = {
    "TsNHNH2": reagent("TsNHNH2", "TsNHNH2", "p-toluenesulfonyl hydrazide"),
    "I2": reagent("I2", "I2", "iodine", "II"), "K2CO3": reagent("K2CO3", "K2CO3", "potassium carbonate", "[K+].[K+].[O-]C([O-])=O"),
    "EtOH": reagent("EtOH", "EtOH", "ethanol", "CCO"), "Ph3PMeBr": reagent("Ph3PMeBr", "Ph3PMeBr", "methyltriphenylphosphonium bromide"),
    "tBuOK": reagent("tBuOK", "tBuOK", "potassium tert-butoxide", "CC(C)(C)[O-].[K+]"),
    "THF": reagent("THF", "THF", "tetrahydrofuran", "C1CCOC1"), "NaBH4": reagent("NaBH4", "NaBH4", "sodium borohydride", "[Na+].[BH4-]"),
    "CeCl3_7H2O": reagent("CeCl3_7H2O", "CeCl3·7H2O", "cerium(III) chloride heptahydrate"),
    "MeOH": reagent("MeOH", "MeOH", "methanol", "CO"), "PhSH": reagent("PhSH", "PhSH", "benzenethiol", "Sc1ccccc1"),
    "CAN": reagent("CAN", "CAN", "ceric ammonium nitrate"), "p_anisidine": reagent("p_anisidine", "p-anisidine", "4-methoxyaniline", "COc1ccc(N)cc1"),
}


def deriv(suffix, substrate, product, reagents, solvent, temp, tim, y, labels, title, special=None):
    evidence = (ev(4),)
    parts = [participant("substrate", substrate, "reactant", evidence)]
    for i, mid in enumerate(reagents):
        parts.append(participant("reagent_" + str(i + 1), mid, "reagent", evidence))
    parts += [participant("product", product, "product", evidence), participant("solvent", solvent, "solvent", evidence)]
    r = reaction(suffix, 4, "specific_experiment", ["derivatization"], title, title, parts,
                 [step("reaction", evidence, [p["participant_id"] for p in parts], temp, tim, special or "reaction under scheme conditions")],
                 [out("yield", "yield", [product], evidence, f"{y}%", {"number": y, "unit": "%"}, "main product",
                      scope="Printed product yield; measurement basis is not stated.")], labels=labels,
                 unresolved=["Yield measurement basis, charges, scale, atmosphere, workup and purification are unreported."])
    if special and "microwave" in special:
        r["special_conditions"]["reactor"] = f("microwave", "microwave reactor (details unreported)", evidence)
    return r


deriv("scheme6a_5a_to_6", M5A, M6, [R["TsNHNH2"], R["I2"], R["K2CO3"]], R["EtOH"],
      {"number": 75, "unit": "°C"}, {"number": 3, "unit": "h"}, 77, ["DB-C-N"], "N-tosylhydrazone formation 5a to 6")
deriv("scheme6a_6_to_7", M6, M7, [R["K2CO3"]], MECHN, {"number": 130, "unit": "°C"},
      {"number": 2, "unit": "h"}, 69, ["C-N", "DB-C-C"], "pyrazole formation 6 to 7", "microwave reaction")
deriv("scheme6a_5a_to_8", M5A, M8, [R["NaBH4"], R["CeCl3_7H2O"]], R["MeOH"], "rt",
      {"number": 2, "unit": "h"}, 93, ["C-H", "O-H"], "chemoselective carbonyl reduction 5a to 8")
r9 = deriv("scheme6a_5a_to_9", M5A, M9, [R["Ph3PMeBr"], R["tBuOK"]], R["THF"], "0 °C to rt",
           {"number": 5, "unit": "h"}, 47, ["DB-C-C"], "Wittig olefination 5a to 9",
           "reaction started at 0 °C and warmed to room temperature; ramp duration not separately reported")
r9["steps"] = [step("cold_start_and_warm", (ev(4),), [p["participant_id"] for p in r9["participants"]],
                    "cold-to-room-temperature program", {"number": 5, "unit": "h"},
                    "combined cold-start/warming stage; the source gives only a total 5 h")]
r9["steps"][0]["temperature"] = f("0 °C to rt", "cold-to-room-temperature program", (ev(4),))
r9["review"]["unresolved"].append("The source does not allocate the 5 h between 0 °C holding and warming/room-temperature holding.")
deriv("scheme6a_5a_to_10", M5A, M10, [R["PhSH"], R["CAN"]], MECHN, {"number": 80, "unit": "°C"},
      {"number": 6, "unit": "h"}, 74, ["C-S", "C-H"], "radical 1,4-addition 5a to 10")
deriv("scheme6b_5j_to_11", molecules["5j"]["molecule_id"], M11, [R["p_anisidine"]], R["MeOH"],
      {"number": 80, "unit": "°C"}, {"number": 1, "unit": "h"}, 87, ["C-N"], "pyrrole formation 5j to 11", "microwave reaction")
deriv("scheme6c_5w_to_12", molecules["5w"]["molecule_id"], M12, [M1], MECHN, "rt", {"number": 2, "unit": "h"},
      63, ["C-C", "C-C"], "cyclobutene formation 5w to 12")
deriv("scheme6c_5w_to_13", molecules["5w"]["molecule_id"], M13, [M1], MECHN, {"number": 80, "unit": "°C"},
      {"number": 2, "unit": "h"}, 57, ["C-C", "C-O"], "chromene formation 5w to 13")
deriv("scheme6d_5x_to_14", molecules["5x"]["molecule_id"], M14, [M1], MECHN, "rt", {"number": 5, "unit": "h"},
      92, ["C-N", "C-C"], "intermolecular heterocyclization 5x to 14")

# Scheme 7 isotope-labeling and atmosphere controls.
H218O = mol("H218O", "H2(18)O", "oxygen-18 water", (ev(4),), "[18OH2]", chirality="achiral",
            geometry="not_applicable", representation="not_applicable")
M18O5A = mol("18O_5a", "18O-5a", "carbonyl-18O-labeled 5a", (ev(4),), chirality="achiral",
             geometry="specified", representation="unresolved",
             stereo_detail={"double_bonds": [{"alkene": "enone C=C", "configuration": "E"}],
                            "assignment": "The source specifies 18O at the carbonyl oxygen."})
parts = pset(M4A, M18O5A, (ev(4),))
parts.insert(2, participant("H218O", H218O, "additive", (ev(4),), equiv=2))
reaction("scheme7a_H218O_labeling", 4, "specific_experiment", ["mechanistic_study"], "Scheme 7a H2(18)O labeling",
         "oxygen-isotope labeling during cascade", parts,
         [step("reaction", (ev(4),), [p["participant_id"] for p in parts], "rt", None,
               "reaction in anhydrous acetonitrile with H2(18)O (2 equiv)")],
         [out("yield", "yield", [M18O5A], (ev(4),), "76%", {"number": 76, "unit": "%"}, "main product",
              scope="Yield of carbonyl-18O-labeled 5a; basis unreported."),
          out("labeling", "other", [M18O5A], (ev(4),), "carbonyl oxygen specifically replaced by 18O",
              "specific carbonyl-18O incorporation", scope="Qualitative isotope-location outcome; quantitation method unreported.")],
         labels=["C-C", "DB-C-O"], unresolved=["Reaction time, charges other than H2(18)O, scale and yield/labeling measurement methods are not reported."])
for tag, atmosphere, y in [("Ar", "Ar atmosphere", 71), ("O2", "O2 balloon", 73)]:
    evidence = (ev(4),)
    parts = pset(M4A, M5A, evidence)
    reaction("scheme7_" + tag + "_control", 4, "specific_experiment", ["mechanistic_study", "control"],
             "Scheme 7 " + atmosphere + " control", "atmosphere-control cascade", parts,
             [step("reaction", evidence, [p["participant_id"] for p in parts], "rt", None,
                   "reaction in acetonitrile", atmosphere=atmosphere)],
             [out("yield", "yield", [M5A], evidence, f"{y}%", {"number": y, "unit": "%"}, "main product",
                  scope="Printed yield; measurement basis unreported.")], labels=["C-C", "DB-C-O"],
             unresolved=["Reaction time, charges, scale and yield measurement basis are not reported.",
                         "The statement that oxygen did not affect the reaction is qualitative; no statistical comparison is provided."])

# Page 5 proposed mechanism (not an experiment) plus intermediates.
M2FPY = mol("2FPy", "2-FPy", "2-fluoropyridine", (ev(5),))
ints = [mol("INT" + str(i), "INT-" + str(i), "proposed intermediate INT-" + str(i), (ev(5),)) for i in range(1, 5)]
parts = [participant("allenols_4", generic_allenol, "reactant", (ev(5),)), participant("reagent_1", M1, "reagent", (ev(5),)),
         participant("water", MWATER, "reagent", (ev(5),)), participant("2FPy", M2FPY, "product", (ev(5),)),
         participant("INT1", ints[0], "unknown", (ev(5),)), participant("INT2", ints[1], "unknown", (ev(5),)),
         participant("INT3", ints[2], "unknown", (ev(5),)), participant("INT4", ints[3], "unknown", (ev(5),)),
         participant("enones_5", generic_enone, "product", (ev(5),))]
reaction("scheme8_proposed_mechanism", 5, "proposed_mechanism", [], "Scheme 8 tentative pathway",
         "proposed electrophilic attack, water addition, proton transfer and dehydration", parts,
         [], [], labels=["C-C", "DB-C-O"],
         unresolved=["This is an author-proposed pathway, not a standalone preparative experiment.",
                     "Intermediate structures and elementary proton/hydrogen transfers are not fully atom-mapped."])

# Source context not represented as preparative reactions.
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:title", "document_id": DOC, "kind": "bibliographic_metadata",
     "evidence_ids": [ev(1)], "value": {"title": "Trifluorosulfonylation Cascade in Allenols: Stereocontrolled Synthesis of Bis(triflyl)enones",
                                                  "doi": "10.1002/chem.202001236"}},
    {"context_id": f"{DOC}:context:dft", "document_id": DOC, "kind": "mechanistic_computation",
     "evidence_ids": [ev(4), ev(5)], "value": {"method": "PCM(MeCN)-M02-6X/6-31+G(d)",
        "reported_summary": "TS-3b 22.7 kcal mol-1; TS-4b 23.3 kcal mol-1; computed C4/C3 kinetic ratio 87:13."}},
    {"context_id": f"{DOC}:context:page6", "document_id": DOC, "kind": "conclusion_and_references",
     "evidence_ids": [ev(6)], "value": "Page 6 contains conclusions, acknowledgements, conflict statement, keywords and bibliography; no new current-paper reaction arrow/table was found."},
    {"context_id": f"{DOC}:context:page7", "document_id": DOC, "kind": "bibliography_and_publication_history",
     "evidence_ids": [ev(7), ev(7, False)],
     "value": {"regions": ["continuation of reference 16", "reference 17", "manuscript publication dates"],
               "finding": "The page contains cited-literature continuation, a cited conformational-preference statement for Tf2C=CH2, and manuscript dates; it contains no current-paper reaction arrow, table, preparative conditions, outcome, or standalone characterization data."}},
    {"context_id": f"{DOC}:context:no_characterization_tables", "document_id": DOC, "kind": "coverage_note",
     "evidence_ids": [ev(i) for i in range(1, 7)], "value": "No standalone MS peak, molecular-weight, fragment-abundance, NMR peak-list or characterization table occurs on pages 1-6."},
]

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/model_extraction.json").write_text(json.dumps({
    "document_id": DOC,
    "pages_attempted": [1, 2, 3, 4, 5, 6],
    "reaction_ids": [r["reaction_id"] for r in reactions],
    "molecule_ids": [m["molecule_id"] for m in bundle["molecules"]],
    "notes": ["Raw machine construction inventory; formal import is runs/current/submitted.json.",
              "No prior workflow or experiment records were used as source data."],
}, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"reactions": len(reactions), "molecules": len(molecules)}, ensure_ascii=False))
