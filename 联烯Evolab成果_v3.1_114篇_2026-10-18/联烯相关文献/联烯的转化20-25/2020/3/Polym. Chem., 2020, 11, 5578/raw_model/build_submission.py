#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for the four-page cycloallene paper."""

from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_625307424a4aa39e"
P = {
    n: {
        "text": f"{DOC}:page{n}:text",
        "image": f"{DOC}:page{n}:image",
    }
    for n in range(1, 5)
}


def f(raw=None, normalized=None, evidence=(), method="original_extraction", status=None,
      derivation=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def direct(raw, normalized, evidence):
    return f(raw, normalized, evidence)


def derived(normalized, evidence, derivation, raw=None):
    return f(raw, normalized, evidence, method="structure_derivation", derivation=derivation)


def inferred(normalized, evidence, derivation, raw=None):
    return f(raw, normalized, evidence, method="inference", derivation=derivation)


def unresolved_raw(raw, evidence):
    return f(raw, None, evidence, status="unreadable")


def set_name(m, original, english, evidence, aliases=None, abbreviations=None):
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = direct(original, original, evidence)
    m["occurrences"][0]["evidence_ids"] = list(evidence)
    m["names"]["original"] = direct(original, original, evidence)
    m["names"]["english"] = direct(english, english, evidence)
    if aliases:
        m["names"]["aliases"] = direct(aliases, aliases, evidence)
    if abbreviations:
        m["names"]["abbreviations"] = direct(abbreviations, abbreviations, evidence)
    m["cas_numbers"] = f()


def concrete_molecule(mid, label, name, smiles, evidence, *, chirality, geometry,
                      representation, stereo_note, aliases=None):
    m = st.new_molecule(mid, DOC, evidence)
    set_name(m, label, name, evidence, aliases=aliases)
    m["representation_type"] = "smiles"
    m["smiles"] = derived(
        smiles, evidence,
        "Candidate connectivity transcribed from the source drawing/name and checked for graph parsing; "
        "parsing does not certify source identity.",
        raw={"source": "source drawing or explicit chemical name", "candidate": smiles},
    )
    diag = st.stereo_inspect(smiles)
    m["stereochemistry"]["chirality_status"] = inferred(
        chirality, evidence, stereo_note,
        raw="Source stereochemical labels and the candidate graph were reviewed."
    )
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, evidence, stereo_note,
        raw="Source stereochemical labels and the candidate graph were reviewed."
    )
    m["stereochemistry"]["representation_status"] = inferred(
        representation, evidence, stereo_note,
        raw="Candidate SMILES was compared with the source drawing."
    )
    m["stereochemistry"]["isomeric_smiles"] = derived(
        diag["canonical_isomeric_smiles"], evidence,
        "RDKit canonical isomeric SMILES for the candidate graph; axial/cyclic stereo limitations remain.",
    )
    m["stereochemistry"]["identity_key"] = derived(
        diag["stereo_graph_key"], evidence,
        "RDKit stereo graph key for comparing this candidate only; it is not a sample identity key.",
    )
    m["identity_note"] = direct(stereo_note, stereo_note, evidence)
    return m, diag


def unresolved_molecule(mid, label, name, evidence, note):
    m = st.new_molecule(mid, DOC, evidence)
    set_name(m, label, name, evidence)
    m["identity_note"] = direct(note, note, evidence)
    for key in ("chirality_status", "geometry_status", "representation_status"):
        m["stereochemistry"][key] = inferred(
            "unresolved", evidence,
            "The source identity/structure is not represented safely enough for a stereochemical claim.",
            raw="No stereochemical claim encoded."
        )
    return m


def template_fact(raw, normalized, evidence):
    return direct(raw, normalized, evidence)


def polymer_molecule(mid, label, name, evidence, description, variables):
    m = st.new_molecule(mid, DOC, evidence)
    set_name(m, label, name, evidence)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": template_fact(
            "Polymer repeat drawing in Scheme 1 or source-defined polymer sample",
            {"description": description}, evidence,
        ),
        "variables": [
            {
                "label": var,
                "attachment_points": template_fact(
                    f"{var} in the source polymer notation", "source-defined polymer variable", evidence),
                "definition": template_fact(
                    f"{var} is not numerically resolved in the supplied pages",
                    "polymerization degree/comonomer count or distribution", evidence),
                "allowed_values": f(),
                "constraints": f(),
            }
            for var in variables
        ],
        "instance_mapping": [],
        "reason": template_fact(
            "Polymeric material with source-defined repeat/distribution",
            "A single finite concrete SMILES would invent chain length, end groups, sequence, or microstructure.",
            evidence,
        ),
    }
    stereo_reason = (
        "The source polymer drawings do not assign chain microstructure, stereocenters, alkene geometry, "
        "or end-group stereochemistry; absence of labels is not treated as achiral or racemic."
    )
    for key, value in (("chirality_status", "unspecified"),
                       ("geometry_status", "unspecified"),
                       ("representation_status", "not_encoded")):
        m["stereochemistry"][key] = inferred(
            value, evidence, stereo_reason, raw="No polymer stereochemical assignment is printed."
        )
    m["identity_note"] = direct(description, description, evidence)
    return m


def participant(pid, mid, role, evidence):
    return st.new_participant(pid, mid, role, evidence)


def base_reaction(rid, evidence, record_type, categories, location, reaction_name,
                  transformation, unresolved):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = direct(categories, categories, evidence)
    r["location"] = direct(location, location, evidence)
    r["reaction_name"] = direct(reaction_name, reaction_name, evidence)
    r["transformation_type"] = inferred(
        transformation, evidence,
        "Transformation class assigned from the source arrow/text; this is not an author-supplied ontology label.",
        raw=reaction_name,
    )
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None,
        "checked_fields": [], "unresolved": list(unresolved),
    }
    r["completeness"] = "partial"
    r["classification"]["limitations"] = [
        "Polymer repeat/end-group atom mapping is unavailable from the supplied pages.",
        "The partial C-C label represents interunit connectivity; chain ends, allene bond-order changes, "
        "hydrogen balance, and catalyst-bound intermediates are unresolved.",
    ]
    return r


def c_c_classification(r, evidence):
    change = {
        "change_id": "b1",
        "atoms": [
            {"atom_map": None, "element": "C", "source_role": "one unsaturated monomer unit"},
            {"atom_map": None, "element": "C", "source_role": "adjacent monomer/growing chain unit"},
        ],
        "before_order": 0,
        "after_order": 1,
        "change_type": "formed",
    }
    why = (
        "Scheme 1 depicts vinyl-addition incorporation into a carbon-chain polymer. Exact repeat-unit atom "
        "correspondence is not supplied, so only the interunit C-C connection is asserted."
    )
    r["classification"]["bond_changes"] = [derived(change, evidence, why)]
    r["classification"]["labels"] = derived(
        [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]}],
        evidence, why,
    )
    r["classification"]["coverage"] = inferred("partial", evidence, why)


def set_step(r, step_id, participants, evidence, operation, temperature=None, time=None,
             atmosphere=None):
    s = st.new_step(step_id, evidence)
    s["participant_ids"] = list(participants)
    s["operation"] = direct(operation, operation, evidence)
    if temperature is not None:
        raw, normalized = temperature
        s["temperature"] = direct(raw, normalized, evidence)
    if time is not None:
        raw, normalized = time
        s["time"] = direct(raw, normalized, evidence)
    if atmosphere is not None:
        s["atmosphere"] = direct(atmosphere, atmosphere, evidence)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, evidence, raw, normalized, *, product_role=None,
            yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = direct(raw, normalized, evidence)
    if product_role is not None:
        o["product_role"] = direct(product_role, product_role, evidence)
    if yield_type is not None:
        o["yield_type"] = direct(yield_type, yield_type, evidence)
    if basis is not None:
        o["measurement_basis"] = direct(basis, basis, evidence)
    if scope is not None:
        o["scope"] = direct(scope, scope, evidence)
    r["outcomes"].append(o)
    return o


def note(r, text, evidence):
    r["notes"].append(direct(text, text, evidence))


dataset = json.loads((RUN / "dataset.json").read_text())
bundle = json.loads((RUN / "submitted.json").read_text())
assert bundle["documents"] == dataset["documents"]

# Molecules and graph diagnostics.
molecules = []
stereo_checks = {}

def add_concrete(*args, **kwargs):
    m, diag = concrete_molecule(*args, **kwargs)
    molecules.append(m)
    stereo_checks[m["molecule_id"]] = diag


add_concrete(
    f"{DOC}:mol:1", "1", "1,2-cyclononadiene", "C1=C=CCCCCCC1", [P[1]["image"], P[2]["text"]],
    chirality="unspecified", geometry="not_applicable", representation="partial",
    stereo_note=("The nine-membered cyclic allene can possess an allene axis, but the paper assigns no axial "
                 "configuration or enantiomeric composition. RDKit does not encode that source-unassigned axis."),
)
add_concrete(
    f"{DOC}:mol:7", "7", "1,2-nonadiene", "C=C=CCCCCCC", [P[1]["image"], P[2]["text"]],
    chirality="achiral", geometry="unspecified", representation="partial",
    stereo_note=("The terminal CH2 end prevents allene axial chirality for this candidate connectivity; the "
                 "source gives no separate alkene-geometry assignment, and toolkit cumulene stereo is limited."),
)
add_concrete(
    f"{DOC}:mol:8", "8", "1,2-cyclotridecadiene", "C1=C=CCCCCCCCCCC1", [P[3]["image"], P[3]["text"]],
    chirality="unspecified", geometry="not_applicable", representation="partial",
    stereo_note=("The thirteen-membered cyclic allene can possess an allene axis, but the paper assigns no axial "
                 "configuration or enantiomeric composition. RDKit does not encode that source-unassigned axis."),
)
add_concrete(
    f"{DOC}:mol:ethylene", "ethylene", "ethylene", "C=C", [P[1]["image"]],
    chirality="achiral", geometry="not_applicable", representation="encoded",
    stereo_note="The source drawing and graph have no stereogenic center or applicable E/Z assignment.",
)
add_concrete(
    f"{DOC}:mol:norbornene", "norbornene drawing", "norbornene", "C1=CC2CCC1C2", [P[1]["image"]],
    chirality="unspecified", geometry="not_applicable", representation="partial",
    stereo_note=("The bicyclic drawing is transcribed as a connectivity candidate; bridgehead stereochemical "
                 "representation is not certified from this compact prior-work graphic."),
)
add_concrete(
    f"{DOC}:mol:cyclooctene", "cyclooctene drawing", "cyclooctene", "C1=CCCCCCC1", [P[1]["image"]],
    chirality="unspecified", geometry="unspecified", representation="partial",
    stereo_note="The compact prior-work drawing does not label cis/trans alkene geometry or enantiomeric form.",
)
add_concrete(
    f"{DOC}:mol:cyclopentene", "cyclopentene drawing", "cyclopentene", "C1=CCCC1", [P[1]["image"]],
    chirality="achiral", geometry="not_applicable", representation="encoded",
    stereo_note="The source drawing and candidate graph contain no assigned or potential independent stereoisomerism.",
)

for mid, label, name, ev, desc, vars_ in [
    ("poly1", "2", "poly(1,2-cyclononadiene)", [P[1]["image"], P[2]["text"]],
     "Vinyl-addition polymer of cyclic allene 1 with an unsaturated ring motif, as drawn in Scheme 1.", ["n"]),
    ("poly7", "poly(7)", "poly(1,2-nonadiene)", [P[1]["image"], P[3]["image"]],
     "Vinyl-addition homopolymer from acyclic allene 7; chain length/end groups are not supplied.", ["n"]),
    ("copoly_generic", "generic poly(1-co-7)", "generic statistical copolymer of 1 and 7", [P[1]["image"]],
     "Generic Scheme 1 copolymer containing cycloallene-derived and acyclic-allene-derived units; no specific feed or product composition is assigned.", ["a", "b", "n"]),
    ("copoly70", "70:30 (1:7) polymer", "statistical copolymer from 70:30 feed of 1 and 7", [P[1]["image"], P[3]["image"]],
     "Statistical copolymer drawn with cycloallene-derived and acyclic-allene-derived units; 70:30 is feed ratio, not assumed measured composition.", ["a", "b", "n"]),
    ("copoly50", "50:50 (1:7) polymer", "statistical copolymer from 50:50 feed of 1 and 7", [P[1]["image"], P[3]["image"]],
     "Statistical copolymer drawn with cycloallene-derived and acyclic-allene-derived units; 50:50 is feed ratio, not assumed measured composition.", ["a", "b", "n"]),
    ("copoly30", "30:70 (1:7) polymer", "statistical copolymer from 30:70 feed of 1 and 7", [P[1]["image"], P[3]["image"]],
     "Statistical copolymer drawn with cycloallene-derived and acyclic-allene-derived units; 30:70 is feed ratio, not assumed measured composition.", ["a", "b", "n"]),
    ("poly8", "poly(8)", "poly(1,2-cyclotridecadiene)", [P[3]["text"], P[3]["image"]],
     "Vinyl-addition polymer from cyclic allene 8; repeat-unit atom mapping, chain length, and end groups are unavailable.", ["n"]),
    ("prev_nbe", "previous-work norbornene/ethylene polymer", "norbornene-ethylene copolymer", [P[1]["image"]],
     "Prior-work cycloolefin/ethylene copolymer with source variables a, b, and n.", ["a", "b", "n"]),
    ("prev_c8", "previous-work cyclooctene/ethylene polymer", "cyclooctene-ethylene copolymer", [P[1]["image"]],
     "Prior-work cyclooctene/ethylene copolymer; sequence and chain length are not supplied.", ["n"]),
    ("prev_c5", "previous-work cyclopentene/ethylene polymer", "cyclopentene-ethylene copolymer", [P[1]["image"]],
     "Prior-work cyclopentene/ethylene copolymer; sequence and chain length are not supplied.", ["n"]),
]:
    molecules.append(polymer_molecule(f"{DOC}:mol:{mid}", label, name, ev, desc, vars_))

for mid, label, name, note_text in [
    ("cat_generic", "Cat", "generic polymerization catalyst", "Scheme 1 uses Cat for the current-work homopolymerization and displays catalysts 3-6 below."),
    ("cat3", "3", "[(pi-allyl)NiOCOCF3]2", "Catalyst 3 is named in the text and drawn in Scheme 1; a safe complete organometallic SMILES was not asserted."),
    ("cat4", "4", "catalyst 4", "The organonickel structure is drawn in Scheme 1, but a complete source-certified SMILES was not asserted."),
    ("cat5", "5", "catalyst 5", "The palladium complex is drawn in Scheme 1, but a complete source-certified SMILES was not asserted."),
    ("cat6", "6", "catalyst 6", "The nickel complex is drawn in Scheme 1, but a complete source-certified SMILES was not asserted."),
    ("agsbf6", "AgSbF6", "silver hexafluoroantimonate", "The activator is explicitly named; no CAS or source structure is supplied."),
    ("M1", "M1", "Fink 1998 catalyst M1", "Prior-work catalyst structure shown only in the compact Scheme 1 inset."),
    ("M2", "M2", "Waymouth 2004 catalyst M2", "Prior-work catalyst structure shown only in the compact Scheme 1 inset."),
    ("MAO", "MAO", "methylaluminoxane", "Named cocatalyst mixture/material; composition is not defined in the supplied pages."),
    ("MMAO", "MMAO", "modified methylaluminoxane", "Named cocatalyst mixture/material; composition is not defined in the supplied pages."),
]:
    molecules.append(unresolved_molecule(
        f"{DOC}:mol:{mid}", label, name,
        [P[1]["image"]] if mid not in ("agsbf6",) else [P[2]["text"], P[2]["image"]], note_text,
    ))

MID = {m["molecule_id"].split(":mol:", 1)[1]: m["molecule_id"] for m in molecules}

reactions = []

# Scheme 1 prior-work arrows.
prior_rows = [
    ("scheme1_prior_norbornene", "norbornene", "prev_nbe", "M1", "MAO"),
    ("scheme1_prior_cyclooctene", "cyclooctene", "prev_c8", "M2", "MMAO"),
    ("scheme1_prior_cyclopentene", "cyclopentene", "prev_c5", "M2", "MMAO"),
]
for rid, mono, product, catalyst, cocat in prior_rows:
    ev = [P[1]["image"]]
    r = base_reaction(
        rid, ev, "literature_example", [], "Page 1, Scheme 1, Previous Work arrow",
        "previous-work cycloolefin/ethylene copolymerization", "vinyl-addition copolymerization",
        ["Scheme 1 provides no experiment-specific quantities, time, temperature, yield, or polymer microstructure."],
    )
    r["participants"] = [
        participant("cyclic_monomer", MID[mono], "reactant", ev),
        participant("ethylene", MID["ethylene"], "reactant", ev),
        participant("catalyst", MID[catalyst], "catalyst", ev),
        participant("cocatalyst", MID[cocat], "reagent", ev),
        participant("polymer", MID[product], "product", ev),
    ]
    c_c_classification(r, ev)
    reactions.append(r)

# Scheme 1 current-work generic arrows.
generic_homo_id = f"{DOC}:scheme1_current_homopolymerization"
ev = [P[1]["image"]]
r = base_reaction(
    "scheme1_current_homopolymerization", ev, "generic_scheme", [],
    "Page 1, Scheme 1, This Work, upper arrow", "vinyl-addition polymerization of cycloallene 1",
    "vinyl-addition homopolymerization",
    ["Generic arrow: exact catalyst, quantities, temperature, time, and yield vary across specific experiments."],
)
r["participants"] = [
    participant("monomer_1", MID["1"], "reactant", ev),
    participant("cat", MID["cat_generic"], "catalyst", ev),
    participant("polymer_2", MID["poly1"], "product", ev),
]
c_c_classification(r, ev)
reactions.append(r)

generic_copoly_id = f"{DOC}:scheme1_current_copolymerization"
r = base_reaction(
    "scheme1_current_copolymerization", ev, "generic_scheme", [],
    "Page 1, Scheme 1, This Work, lower arrow", "copolymerization of cyclic allene 1 and acyclic allene 7",
    "vinyl-addition copolymerization",
    ["Generic arrow: feed ratio, catalyst amount, temperature, time, and yield are experiment-specific or unreported."],
)
r["participants"] = [
    participant("monomer_1", MID["1"], "reactant", ev),
    participant("monomer_7", MID["7"], "reactant", ev),
    participant("catalyst_3", MID["cat3"], "catalyst", ev),
    participant("copolymer", MID["copoly_generic"], "product", ev),
]
note(r, "The displayed copolymer is generic; it is not assigned the 70:30 sample composition by the source.", ev)
c_c_classification(r, ev)
reactions.append(r)

# Table 1: eight independent optimization rows.
table1 = [
    (1, "cat3", False, 0.1, 1, 23, 24400, 35600, 1.11, 94),
    (2, "cat3", False, 0.1, 1, 35, 24400, 50500, 1.11, 93),
    (3, "cat3", False, 0.1, 1, 50, 24400, 48900, 1.13, 97),
    (4, "cat3", False, 0.1, 1, 65, 24400, 73300, 1.12, 53),
    (5, "cat4", False, 1.5, 1, 23, 24400, 133000, 1.62, 57),
    (6, "cat5", True, 1.5, 24, 23, 24400, 32000, 1.37, 63),
    (7, "cat6", True, 1.5, 1, 23, 24400, None, None, 59),
    (8, "cat3", True, 1.5, 24, 23, 24400, 8100, 2.73, 99),
]
for entry, cat, activated, conc, hours, temp, mtheo, mexpt, dispersity, yield_pct in table1:
    ev = [P[2]["image"], P[2]["text"]]
    unresolved = [
        "The [1]0 column gives a number but the supplied table header does not print a concentration unit.",
        "Monomer/catalyst charges, catalyst loading, solvent, atmosphere, pressure, workup, and purification are unavailable without SI.",
        "Yield measurement basis is not stated in the supplied main text/table.",
        "Polymer chain length, end groups, microstructure, and complete atom mapping remain unresolved.",
    ]
    if entry == 7:
        unresolved.append("The material is reported insoluble; experimental Mn and dispersity are table dashes.")
    r = base_reaction(
        f"table1_entry_{entry}", ev, "specific_experiment", ["optimization"],
        f"Page 2, Table 1, entry {entry}", "vinyl-addition homopolymerization of 1",
        "vinyl-addition homopolymerization", unresolved,
    )
    r["general_scheme_ids"] = [generic_homo_id]
    r["participants"] = [
        participant("monomer_1", MID["1"], "reactant", ev),
        participant("catalyst", MID[cat], "catalyst", ev),
        participant("polymer_2", MID["poly1"], "product", ev),
    ]
    if activated:
        r["participants"].append(participant("AgSbF6", MID["agsbf6"], "additive", ev))
    r["solvent_conditions"]["concentration"] = unresolved_raw(
        {"table_header": "[1]0", "cell": str(conc), "unit": "not printed"}, ev)
    pids = [p["participant_id"] for p in r["participants"]]
    set_step(
        r, "polymerization", pids, ev, "polymerization under Table 1 conditions",
        temperature=(f"{temp} °C", {"number": temp, "unit": "°C"}),
        time=(f"{hours} h", {"number": hours, "unit": "h"}),
    )
    outcome(r, "yield", "yield", [MID["poly1"]], ev, f"{yield_pct}%", {"number": yield_pct, "unit": "%"},
            product_role="polymer product")
    theo_basis = ("calculated using [1]0:[Metal]0 × 122.12" if entry != 5
                  else "calculated using [1]0:[Catalyst]0 × 122.12")
    outcome(r, "Mn_theoretical", "other", [MID["poly1"]], ev, f"{mtheo} Da", f"{mtheo} Da",
            basis=theo_basis, scope="theoretical number-average molecular weight")
    if mexpt is not None:
        outcome(r, "Mn_experimental", "other", [MID["poly1"]], ev, f"{mexpt} Da", f"{mexpt} Da",
                basis="GPC (triple detection)", scope="experimental number-average molecular weight")
        outcome(r, "dispersity", "other", [MID["poly1"]], ev, str(dispersity), dispersity,
                basis="calculated using Mw/Mn", scope="Đ")
    else:
        outcome(r, "insoluble_material", "other", [MID["poly1"]], ev, "Insoluble material", "insoluble material",
                basis="Table 1 footnote e", scope="reason Mn and Đ are not reported")
    if activated:
        note(r, "Table 1 footnote f: activated with AgSbF6.", ev)
    if entry == 1:
        note(r, "The prose identifies this 94% material as polymer 2, a white powder, and reports a recurring high-Mn impurity.", ev)
    c_c_classification(r, ev)
    reactions.append(r)

# Two kinetic benchmark polymerizations, distinct from Table 1 rows.
kinetics = [
    ("kinetics_monomer_1", "1", "poly1", "within 20 min", {"number": 20, "unit": "min", "qualifier": "le"}, 0.296),
    ("kinetics_monomer_7", "7", "poly7", "approximately 60 min", {"number": 60, "unit": "min", "qualifier": "approx"}, 0.052),
]
for rid, mono, product, raw_time, norm_time, rate in kinetics:
    ev = [P[2]["text"], P[2]["image"]]
    r = base_reaction(
        rid, ev, "specific_experiment", ["mechanistic_study"],
        "Page 2, kinetic comparison prose and Fig. 1", f"kinetic homopolymerization of monomer {mono}",
        "vinyl-addition homopolymerization kinetics",
        ["Exact quantities, catalyst loading, concentration unit, temperature, atmosphere, and ESI procedure are unavailable.",
         "The printed kobs values have no unit in the supplied main-text prose/figure."],
    )
    if mono == "1":
        r["general_scheme_ids"] = [generic_homo_id]
    r["participants"] = [
        participant(f"monomer_{mono}", MID[mono], "reactant", ev),
        participant("catalyst_3", MID["cat3"], "catalyst", ev),
        participant("polymer", MID[product], "product", ev),
    ]
    set_step(r, "kinetic_polymerization", [p["participant_id"] for p in r["participants"]], ev,
             "polymerization monitored by periodic 1H NMR spectroscopy", time=(raw_time, norm_time))
    outcome(r, "conversion", "conversion", [MID[mono]], ev, ">95%", {"number": 95, "unit": "%", "qualifier": "gt"},
            basis="1H NMR spectroscopy", scope="monomer consumption/conversion")
    outcome(r, "kobs", "other", [MID[mono]], ev, f"kobs = {rate}", {"number": rate, "unit_reported": False},
            basis="slope of Fig. 1 ln(M0/Mt) versus time plot", scope="observed polymerization rate constant")
    c_c_classification(r, ev)
    reactions.append(r)

# Living-character aliquot experiment.
ev = [P[3]["text"], P[3]["image"]]
r = base_reaction(
    "living_conversion_series", ev, "specific_experiment", ["mechanistic_study"],
    "Page 3, living-character prose and Fig. 2B", "conversion-dependent polymerization study of 1",
    "living polymerization conversion series",
    ["Aliquot times, exact conversion/Mn data points, catalyst loading, temperature, and full ESI procedure are not tabulated in the supplied pages."],
)
r["general_scheme_ids"] = [generic_homo_id]
r["participants"] = [
    participant("monomer_1", MID["1"], "reactant", ev),
    participant("catalyst_3", MID["cat3"], "catalyst", ev),
    participant("polymer_2", MID["poly1"], "product", ev),
]
set_step(r, "polymerization_sampling", [p["participant_id"] for p in r["participants"]], ev,
         "aliquots periodically removed during polymerization and quenched by exposure to air")
outcome(r, "Mn_conversion_relation", "other", [MID["poly1"]], ev,
        "linear correlation between Mn and monomer conversion", "linear correlation",
        basis="conversion by 1H NMR spectroscopy; Mn by GPC", scope="living-character diagnostic")
outcome(r, "dispersity_bound", "other", [MID["poly1"]], ev, "Đ ≤ 1.2",
        {"number": 1.2, "qualifier": "le"},
        basis="GPC", scope="dispersity during conversion series")
c_c_classification(r, ev)
reactions.append(r)

# Chain extension is a separate experiment, not merged with the conversion series.
r = base_reaction(
    "chain_extension_1", ev, "specific_experiment", ["mechanistic_study"],
    "Page 3, Fig. 2A and chain-extension prose", "chain extension via sequential addition of 1",
    "living-polymer chain extension",
    ["Sequential-addition quantities, intervals, temperature, conversion, and detailed procedure are only referenced to unavailable ESI."],
)
r["general_scheme_ids"] = [generic_homo_id]
r["participants"] = [
    participant("monomer_1", MID["1"], "reactant", ev),
    participant("catalyst_3", MID["cat3"], "catalyst", ev),
    participant("chain_extended_polymer", MID["poly1"], "product", ev),
]
set_step(r, "sequential_addition", [p["participant_id"] for p in r["participants"]], ev,
         "sequential addition of monomer 1 to the living polymerization")
outcome(r, "chain_extension", "other", [MID["poly1"]], ev, "chain extension was successful", True,
        basis="GPC trace shift in Fig. 2A", scope="qualitative chain-extension result")
c_c_classification(r, ev)
reactions.append(r)

# Table 2: five independent (co)polymerizations.
table2 = [
    (1, 100, 0, "poly1", 76, 61356, 1.32),
    (2, 0, 100, "poly7", 74, 67570, 1.15),
    (3, 70, 30, "copoly70", 97, 129100, 1.38),
    (4, 50, 50, "copoly50", 92, 98360, 1.85),
    (5, 30, 70, "copoly30", 78, 99270, 1.75),
]
for entry, ratio1, ratio7, product, yield_pct, mn, disp in table2:
    ev = [P[3]["image"], P[3]["text"]]
    r = base_reaction(
        f"table2_entry_{entry}", ev, "specific_experiment", ["substrate_scope"],
        f"Page 3, Table 2, entry {entry}", "(co)polymerization of 1 and/or 7",
        "vinyl-addition (co)polymerization",
        ["Table 2 reports feed ratio but not measured copolymer composition.",
         "Catalyst charge/loading, reaction time, exact temperature, concentration, solvent, atmosphere, workup, and purification are unavailable without ESI.",
         "Yield measurement basis is not stated in the supplied main text/table."],
    )
    if ratio1 == 100:
        r["general_scheme_ids"] = [generic_homo_id]
    elif ratio1 and ratio7:
        r["general_scheme_ids"] = [generic_copoly_id]
    parts = []
    if ratio1:
        parts.append(participant("monomer_1", MID["1"], "reactant", ev))
    if ratio7:
        parts.append(participant("monomer_7", MID["7"], "reactant", ev))
    parts += [participant("catalyst_3", MID["cat3"], "catalyst", ev),
              participant("polymer", MID[product], "product", ev)]
    r["participants"] = parts
    set_step(r, "polymerization", [p["participant_id"] for p in parts], ev,
             f"polymerization with feed ratio 1:7 = {ratio1}:{ratio7}")
    note(r, f"Source feed ratio 1:7 = {ratio1}:{ratio7}; this is not silently converted to product composition.", ev)
    outcome(r, "yield", "yield", [MID[product]], ev, f"{yield_pct}%", {"number": yield_pct, "unit": "%"},
            product_role="polymer product")
    outcome(r, "Mn", "other", [MID[product]], ev, f"{mn} Da", f"{mn} Da",
            basis="GPC", scope="number-average molecular weight")
    outcome(r, "dispersity", "other", [MID[product]], ev, str(disp), disp,
            basis="calculated using Mw/Mn", scope="Đ")
    c_c_classification(r, ev)
    reactions.append(r)

# Larger-ring monomer 8 experiment.
ev = [P[3]["text"], P[3]["image"]]
r = base_reaction(
    "polymerization_monomer_8", ev, "specific_experiment", ["substrate_scope"],
    "Page 3, lower-right prose", "polymerization of 1,2-cyclotridecadiene (8)",
    "vinyl-addition homopolymerization",
    ["The source says optimized polymerization conditions but delegates exact quantities, temperature, solvent, atmosphere, workup, and purification to unavailable ESI.",
     "Conversion measurement method is not reported in the supplied pages.",
     "No isolated yield is reported."],
)
r["participants"] = [
    participant("monomer_8", MID["8"], "reactant", ev),
    participant("catalyst_3", MID["cat3"], "catalyst", ev),
    participant("polymer_8", MID["poly8"], "product", ev),
]
set_step(r, "polymerization", [p["participant_id"] for p in r["participants"]], ev,
         "polymerization under the authors' optimized conditions (exact ESI conditions unavailable)",
         time=("72 h", {"number": 72, "unit": "h"}))
outcome(r, "conversion", "conversion", [MID["8"]], ev, "67%", {"number": 67, "unit": "%"},
        scope="monomer conversion")
outcome(r, "Mn", "other", [MID["poly8"]], ev, "Mn = 47 289", "47289 Da",
        basis="high-temperature GPC", scope="number-average molecular weight")
outcome(r, "dispersity", "other", [MID["poly8"]], ev, "Đ = 1.99", 1.99,
        basis="high-temperature GPC", scope="Đ")
note(r, "The product was poorly soluble in CHCl3, DCM, THF, PhMe, etc.; a portion dissolved in CDCl3 and its 1H NMR spectrum was consistent with the proposed polymer.", ev)
c_c_classification(r, ev)
reactions.append(r)

bundle["molecules"] = molecules
bundle["reactions"] = reactions
bundle["source_context"] = [
    {
        "context_id": f"{DOC}:context:article_metadata",
        "document_id": DOC,
        "kind": "bibliographic_metadata",
        "evidence_ids": [P[1]["text"], P[1]["image"]],
        "title": "Vinyl-addition polymerizations of cycloallenes: synthetic access to congeners of cyclic-olefin polymers",
        "authors": ["Nicholas J. Galan", "Justin M. Burroughs", "Christopher R. Maroon", "Brian K. Long", "Johnathan N. Brantley"],
        "journal": "Polymer Chemistry", "year": 2020, "volume": 11, "pages": "5578-5581",
        "doi": "10.1039/d0py01018a",
    },
    {
        "context_id": f"{DOC}:context:polymer2_characterization",
        "document_id": DOC, "kind": "characterization_and_material_property",
        "evidence_ids": [P[2]["text"], P[2]["image"]],
        "observations": [
            "Table 1 entry 1 material: Mn 35.6 kDa and Đ 1.11 by GPC; a recurring high-Mn impurity was observed.",
            "Thermogravimetric analysis gave Td = 350 °C.",
            "No thermal transitions were observed by differential scanning calorimetry.",
            "Predicted ring strain of monomer 1 is approximately 12 kcal mol-1.",
        ],
        "not_a_separate_synthesis": True,
    },
    {
        "context_id": f"{DOC}:context:tensile_testing",
        "document_id": DOC, "kind": "film_preparation_and_mechanical_testing",
        "evidence_ids": [P[3]["text"], P[3]["image"]],
        "observations": [
            "Materials were cast into films from CHCl3 for tensile testing.",
            "70:30 feed material: Young's modulus 1367 ± 15 MPa; strain at break 2.5% ± 0.6%.",
            "50:50 feed material: Young's modulus 263 ± 30 MPa; strain at break 5.2% ± 2.4%.",
            "Homopolymer 1 was too brittle; homopolymer 7 and the 30 mol% 1 copolymer were too tacky for mechanical analysis.",
        ],
        "not_a_synthetic_reaction": True,
    },
    {
        "context_id": f"{DOC}:context:footnote24_unavailable_si",
        "document_id": DOC, "kind": "unresolved_mechanistic_troubleshooting",
        "evidence_ids": [P[4]["text"], P[4]["image"]],
        "observations": [
            "Footnote 24 says a series of experiments tested oxidative chain-coupling or catalyst chain transfer as origins of the high-Mn species.",
            "The studies did not identify the origin; experiment count and conditions are only in unprovided ESI and are not fabricated here.",
        ],
        "record_created": False,
    },
    {
        "context_id": f"{DOC}:context:si_status",
        "document_id": DOC, "kind": "source_limitation",
        "evidence_ids": [P[1]["text"], P[1]["image"]],
        "observation": "The article advertises ESI, but no SI file was supplied or explicitly associated in this workspace.",
    },
]

assert len(reactions) == 23
assert sum(r["record_type"] == "specific_experiment" for r in reactions) == 18
st.write(RUN / "submitted.json", bundle)
st.write(RUN / "raw_model/stereo_checks.json", stereo_checks)

mentions = st.source_mentions(dataset, ["1,2-cyclononadiene", "1,2-nonadiene", "1,2-cyclotridecadiene", "AgSbF6"])
st.write(RUN / "raw_model/source_mentions.json", mentions)

inventory = {
    "document_id": DOC,
    "pages_attempted": [1, 2, 3, 4],
    "reaction_records": len(reactions),
    "specific_experiments": 18,
    "generic_current_work_schemes": 2,
    "literature_examples": 3,
    "molecule_records": len(molecules),
    "page_inventory": {
        "1": {"scheme_arrows": 5, "records": [r["reaction_id"] for r in reactions[:5]]},
        "2": {"table1_rows": 8, "kinetic_experiments": 2,
              "records": [r["reaction_id"] for r in reactions[5:15]]},
        "3": {"living_conversion_series": 1, "chain_extension": 1, "table2_rows": 5, "large_ring_experiment": 1,
              "records": [r["reaction_id"] for r in reactions[15:23]]},
        "4": {"new_synthesis_records": 0, "context": [f"{DOC}:context:footnote24_unavailable_si"]},
    },
    "semantic_accuracy_verified": False,
    "chemical_accuracy_verified": False,
}
st.write(RUN / "raw_model/pages_0001_0004_extraction.json", inventory)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1, "state": "attempted",
            "regions": [
                {"region": "Scheme 1, Previous Work panel", "finding": "Three cited cycloolefin/ethylene polymerization arrows recorded as literature examples, not this paper's experiments.", "evidence_ids": [P[1]["image"]]},
                {"region": "Scheme 1, This Work panel", "finding": "Two current-work generic arrows: homopolymerization of 1 and copolymerization of 1 with 7; catalyst structures 3-6 inventoried.", "evidence_ids": [P[1]["image"]]},
                {"region": "Title, abstract, introduction, and ESI footnote", "finding": "Bibliographic title/DOI captured; ESI is advertised but no associated SI file was provided.", "evidence_ids": [P[1]["text"], P[1]["image"]]},
            ],
            "record_ids": [r["reaction_id"] for r in reactions[:5]],
            "unresolved": ["Compact organometallic catalyst drawings were not converted to source-certified complete SMILES.", "Generic polymer chain lengths/end groups and stereochemistry are unreported."],
        },
        {
            "page": 2, "state": "attempted",
            "regions": [
                {"region": "Main-text model-substrate and initial-result prose", "finding": "Table 1 entry 1 prose merged as evidence for the same 94% experiment; characterization retained as context.", "evidence_ids": [P[2]["text"], P[2]["image"]]},
                {"region": "Fig. 1 and kinetic prose", "finding": "Two independent kinetic polymerizations (1 and 7) recorded with >95% conversion times, kobs values, and 1H NMR basis.", "evidence_ids": [P[2]["text"], P[2]["image"]]},
                {"region": "Table 1, entries 1-8 and footnotes a-g", "finding": "Eight condition-specific optimization records; AgSbF6 added for footnote-f rows; insoluble entry 7 retains table dashes for Mn/Đ.", "evidence_ids": [P[2]["image"], P[2]["text"]]},
            ],
            "record_ids": [r["reaction_id"] for r in reactions[5:15]],
            "unresolved": ["[1]0 unit is absent from the printed header.", "Yield measurement basis, exact charges, solvent, atmosphere, workup, and purification are unavailable without SI.", "kobs units are not printed."],
        },
        {
            "page": 3, "state": "attempted",
            "regions": [
                {"region": "Living-character prose and Fig. 2", "finding": "Aliquot/conversion study and sequential-addition chain extension recorded as two distinct mechanistic-study experiments.", "evidence_ids": [P[3]["text"], P[3]["image"]]},
                {"region": "Table 2, entries 1-5", "finding": "Five independent homo/copolymerization records with feed ratios, yields, Mn, and Đ; feed ratio not reinterpreted as measured composition.", "evidence_ids": [P[3]["image"], P[3]["text"]]},
                {"region": "Fig. 3 and tensile-testing prose", "finding": "Film casting and tensile results classified as material testing in source_context, not synthetic reactions.", "evidence_ids": [P[3]["text"], P[3]["image"]]},
                {"region": "Lower-right monomer 8 prose", "finding": "One 72 h polymerization of 8 recorded with 67% conversion, high-temperature GPC Mn/Đ, and solubility/NMR notes.", "evidence_ids": [P[3]["text"], P[3]["image"]]},
            ],
            "record_ids": [r["reaction_id"] for r in reactions[15:23]],
            "unresolved": ["Table 2 exact catalyst loading, time, temperature, concentration, solvent, and yield basis are not supplied.", "Fig. 2 numeric series and chain-extension charges are only detailed in unavailable ESI.", "Conversion method for monomer 8 is unreported."],
        },
        {
            "page": 4, "state": "attempted",
            "regions": [
                {"region": "Conclusion", "finding": "General current-work interpretation only; no new experiment-specific arrow or table row.", "evidence_ids": [P[4]["text"], P[4]["image"]]},
                {"region": "Footnote/reference 24", "finding": "Mentions an unresolved series of high-Mn impurity mechanism tests in ESI; insufficient detail/count to create experimental records.", "evidence_ids": [P[4]["text"], P[4]["image"]]},
                {"region": "References 1-31", "finding": "Cited literature list only; not imported as experiments newly reported in this paper.", "evidence_ids": [P[4]["image"]]},
            ],
            "record_ids": [],
            "unresolved": ["Footnote 24 experiment identities, conditions, and number are unavailable because SI was not provided."],
        },
    ],
}
st.write(RUN / "coverage/pages_0001_0004.json", coverage)

tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    if task["document_id"] == DOC and task["task_id"] in {P[n]["text"] for n in range(1, 5)}:
        task["state"] = "attempted"
        task["reason"] = "Original page image and extracted text inspected; detailed regions and limitations are in coverage/pages_0001_0004.json."
st.write(RUN / "tasks.json", tasks)

workflow = f"""# Extraction workflow log

## Batch and source scope

- Batch: `allene-stereo-v31-200-20261017`.
- Document: `{DOC}`.
- Source mapping: `Polym. Chem., 2020, 11, 5578.pdf`.
- DOI: `10.1039/d0py01018a`.
- Title: *Vinyl-addition polymerizations of cycloallenes: synthetic access to congeners of cyclic-olefin polymers*.
- Attempted pages: original PDF pages 1-4, exactly four pages.

`runs/current/source_mapping.json` identifies an original journal PDF (`presentation_conversion: false`). The article advertises ESI, but no SI PDF was supplied or explicitly associated; `supplementary_material.status` remains `unknown`. Exact SI procedures, workups, purification, chain-extension quantities, and footnote-24 mechanism tests were not fetched or inferred.

The pre-existing `runs/workflow.md` described a different paper and was not used as source evidence. This log replaces it with document-specific decisions.

## Instructions and visual review

Read `AGENTS.md`, the complete `skills/chemical-reaction-extraction/SKILL.md`, and the required v3.1 references: `extraction_standard.md`, `stereochemistry.md`, `structured_fields.md`, `quality_rules.md`, `literature.schema.json`, and `bond_taxonomy.json`. No extraction tool, schema, scientific-quality rule, test, or control program was modified.

All four rendered original pages were viewed at original detail. Scheme 1, Table 1, Fig. 1, Table 2, and Fig. 2 were additionally inspected using source-lineage crops saved under `runs/current/raw_model/`, each with a `.evidence.json` receipt. Extracted text was used for prose/table cross-checking and literal source mentions, not as a substitute for viewing page images.

The coverage contract is `runs/current/coverage/pages_0001_0004.json`: exactly pages 1-4, every page `attempted`, every page has evidence-linked regions, and `semantic_accuracy_verified` remains false.

## Page inventory and record decisions

### Page 1

- Scheme 1 contains five arrows. Three arrows in “Previous Work” are separate `literature_example` records; cited work is not counted as this paper's experiments.
- Two arrows in “This Work” are `generic_scheme` records: homopolymerization of cyclic allene 1 and copolymerization of 1 with acyclic allene 7.
- Catalysts 3-6 and prior-work M1/M2 are registered from the source labels/drawings. Dense organometallic structures remain unresolved rather than receiving unsafe complete SMILES.

### Page 2

- Table 1 contributes eight independent optimization experiments. Different catalysts, temperature, time, concentration cells, and AgSbF6 activation are not merged.
- Entry 1 and the adjacent prose's 94% white polymer 2 are the same experiment and are merged as evidence, not duplicated.
- Footnote f is structured by adding AgSbF6 as an additive participant to entries 6-8. Entry 7 retains insolubility and table dashes for Mn/Đ.
- The `[1]0` header supplies no printed unit. The numeric cells are retained as unresolved raw concentration expressions; no M unit is invented.
- Fig. 1/prose supplies two separate kinetic experiments for 1 and 7. Both >95% conversion observations are targeted to their monomers with the explicit 1H NMR method. Printed `kobs` values are retained with unit unreported.
- Yield basis is not explicitly stated in Table 1/main text, so each yield keeps a null measurement basis and an unresolved review item.

### Page 3

- Periodic aliquot/air-quench conversion tracking and sequential-addition chain extension are two distinct mechanistic-study experiments. Their ESI-only quantities and time points remain missing.
- Table 2 contributes five independent (co)polymerizations. The `1:7` values are feed ratios and are not silently rewritten as measured product composition. Mn is GPC-measured; Đ is calculated from Mw/Mn; yield basis is unreported.
- Film casting, tensile testing, Young's modulus, strain at break, brittleness, and tackiness are sample preparation/material characterization, saved to `source_context`, not synthesis records.
- The 72 h polymerization of 8 is one scope experiment with 67% conversion. Conversion method and isolated yield are unreported. High-temperature GPC Mn/Đ and solubility/NMR observations are retained.

### Page 4

- Conclusion and references contain no new specific reaction arrows or table rows.
- Footnote 24 describes a *series* of experiments about the high-Mn impurity but gives neither count nor conditions outside unavailable ESI. It is retained in `source_context`; no fabricated single experiment stands in for an unknown series.
- References were not imported as new experiments.

## Molecules and stereochemistry

- Registry entities: {len(molecules)} molecules/templates/catalysts/reagents.
- Candidate SMILES for 1, 7, 8, ethylene, norbornene, cyclooctene, and cyclopentene were checked with `standard stereo`; diagnostics are in `runs/current/raw_model/stereo_checks.json`.
- Cyclic allenes 1 and 8 retain potential but source-unassigned axial chirality (`chirality_status=unspecified`, `representation_status=partial`). No racemate, enantiomer, or axial assignment is invented.
- Terminal allene 7 is structurally achiral at the allene axis, but cumulene/double-bond representation remains partial because toolkit stereo diagnostics do not certify source geometry.
- Polymer products are generic templates. Chain length, end groups, sequence, alkene geometry, and stereochemical microstructure remain unreported (`chirality_status=unspecified`, `geometry_status=unspecified`, `representation_status=not_encoded`).
- No ee, er, dr, E/Z ratio, absolute configuration, or CAS number is reported; none is invented.
- All imported records remain `needs_review`; source-reviewed count is zero.

## Classification

- Successful polymerization records have partial `C-C` classification for interunit carbon-carbon connectivity.
- Exact repeat-unit atom mapping, allene bond-order changes, end-group hydrogen balance, catalyst-bound atoms, and polymer sequence are unavailable. No classification is marked complete, no unchanged bond is mislabeled, and no model votes are fabricated.

## Counts and deduplication

- Total reaction records: 23.
- Specific experiments: 18 (Table 1: 8; kinetic comparison: 2; living conversion series: 1; chain extension: 1; Table 2: 5; monomer 8: 1).
- Non-specific records: 5 (two current-work generic schemes and three prior-work literature examples).
- Repeated prose/table depictions of Table 1 entry 1 were merged. Different table rows and kinetic/material-property studies were not merged by product label.

## Raw outputs and formal import

- Reproducible builder: `runs/current/raw_model/build_submission.py`.
- Raw page inventory: `runs/current/raw_model/pages_0001_0004_extraction.json`.
- Literal source mentions: `runs/current/raw_model/source_mentions.json`.
- Stereo diagnostics: `runs/current/raw_model/stereo_checks.json`.
- Source-lineage crops/receipts: `runs/current/raw_model/page*.jpg*`.
- Formal machine submission: `runs/current/submitted.json`.

## Known unresolved items

- No SI was supplied: exact charges, catalyst loadings, solvent, atmosphere, workup, purification, chain-extension schedule, and footnote-24 tests remain unavailable.
- Table 1 `[1]0` unit and printed `kobs` units are absent.
- Yield measurement basis is absent for Table 1 and Table 2; it is not assumed to be isolated or NMR yield.
- Table 2 exact reaction times/temperatures and measured copolymer composition are not supplied.
- Monomer-8 conversion method and isolated yield are unreported.
- Organometallic catalyst structures and polymer atom mappings need expert review.
- `semantic_accuracy_verified` and `chemical_accuracy_verified` remain false.

## Validation result

Validation, merge, audit, and test results are appended after the formal tool runs.
"""
(ROOT / "runs/workflow.md").write_text(workflow)

print(json.dumps({"molecules": len(molecules), "reactions": len(reactions), "specific_experiments": 18}, ensure_ascii=False))
