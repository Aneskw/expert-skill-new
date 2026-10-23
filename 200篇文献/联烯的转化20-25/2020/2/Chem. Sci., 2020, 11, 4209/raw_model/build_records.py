#!/usr/bin/env python3
"""Build source-derived v3.1 records for PDF pages 1--6.

This file is a transparent extraction recipe, not a chemistry validator.  Exact
SMILES for the densely fused source drawings are deliberately not asserted.
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
DOC = "doc_14722b08ab785329"


def ev(*pages):
    return [f"{DOC}:page{p}:text" for p in pages]


def F(raw=None, normalized=None, pages=(), method="original_extraction", derivation=None,
      status=None):
    if normalized is None and raw is not None:
        normalized = raw
    return st.fact(raw, normalized, ev(*pages), method=method, derivation=derivation,
                   status=status)


def quantity(number, unit, raw_text=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if raw_text is not None:
        value["raw_text"] = raw_text
    if qualifier is not None:
        value["qualifier"] = qualifier
    return value


def molecule(mid, label, pages, *, name=None, chirality="unresolved",
             geometry="unresolved", representation="not_encoded",
             assignment=None, mixture=None, note=None):
    m = st.new_molecule(mid, DOC, ev(*pages))
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, pages)
    m["names"]["original"] = F(label, label, pages)
    if name:
        m["names"]["english"] = F(name, name, pages)
    m["representation_type"] = "unresolved"
    m["smiles"] = F()
    m["generic_template"] = None
    m["stereochemistry"]["chirality_status"] = F(
        f"source status: {chirality}", chirality, pages, method="inference",
        derivation="Status assigned conservatively from the supplied drawing/text; exact structure is not encoded.")
    m["stereochemistry"]["geometry_status"] = F(
        f"source status: {geometry}", geometry, pages, method="inference",
        derivation="Status assigned conservatively from the supplied drawing/text; exact structure is not encoded.")
    m["stereochemistry"]["representation_status"] = F(
        f"source status: {representation}", representation, pages, method="inference",
        derivation="The source drawing was inspected but no source-verified isomeric SMILES is available.")
    if assignment:
        m["stereochemistry"]["assignment"] = F(assignment, assignment, pages)
    if mixture:
        m["stereochemistry"]["mixture"] = F(mixture, mixture, pages)
    m["identity_note"] = F(
        note or "Identity retained by paper label and source drawing; exact connectivity/stereochemistry was not promoted to SMILES.",
        note or "Identity retained by paper label and source drawing; exact connectivity/stereochemistry was not promoted to SMILES.",
        pages)
    return m


def generic_molecule(mid, label, pages, fixed, variables, note):
    m = molecule(mid, label, pages, note=note)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": F(fixed, fixed, pages),
        "variables": [],
        "instance_mapping": [],
        "reason": F(note, note, pages),
    }
    for variable in variables:
        m["generic_template"]["variables"].append({
            "label": variable,
            "attachment_points": F(f"attachment of {variable} shown in source drawing",
                                     f"attachment of {variable} shown in source drawing", pages),
            "definition": F(f"variable {variable}", f"variable {variable}", pages),
            "allowed_values": F(),
            "constraints": F(),
        })
    return m


def reaction(suffix, pages, *, record_type="specific_experiment", categories=None,
             name=None, transformation=None, location=None, unresolved=None):
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, ev(*pages))
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories or ["other"], categories or ["other"], pages)
    r["location"] = F(location or suffix, {"pages": list(pages), "region": location or suffix}, pages)
    if name:
        r["reaction_name"] = F(name, name, pages)
    if transformation:
        r["transformation_type"] = F(transformation, transformation, pages)
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None,
        "checked_fields": [],
        "unresolved": list(unresolved or [
            "Exact source structures are retained by labels/drawings and need expert structure review.",
            "No explicitly supplied SI was associated; unreported experimental details remain missing.",
        ]),
    }
    r["completeness"] = "partial"
    return r


def participant(pid, mid, role, pages, *, mol_percent=None, equivalents=None,
                reference=None, product_role=None, stage="reaction"):
    p = st.new_participant(pid, mid, role, ev(*pages), stage=stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, pages)
    if mol_percent is not None:
        p["charges"]["mol_percent"] = F(
            f"{mol_percent} mol%", quantity(mol_percent, "mol%", f"{mol_percent} mol%"), pages)
        p["equivalence_reference"] = F(f"relative to {reference}", reference, pages)
    if equivalents is not None:
        p["charges"]["equivalents"] = F(
            f"{equivalents} equiv", quantity(equivalents, "equiv", f"{equivalents} equiv"), pages)
        p["equivalence_reference"] = F(f"relative to {reference}", reference, pages)
    return p


def set_solvent(r, mid, name, pages, concentration=None, basis=None):
    r["participants"].append(participant(f"solvent_{len(r['participants'])+1}", mid, "solvent", pages))
    r["solvent_conditions"]["composition"] = F(name, [name], pages)
    if concentration is not None:
        r["solvent_conditions"]["concentration"] = F(
            f"{concentration} M", quantity(concentration, "M", f"{concentration} M"), pages)
        r["solvent_conditions"]["concentration_basis"] = F(
            f"concentration relative to {basis}", basis, pages)


def add_step(r, sid, pages, *, participant_ids=None, operation=None, temperature=None,
             time=None, atmosphere=None, pressure=None):
    s = st.new_step(sid, ev(*pages))
    s["participant_ids"] = list(participant_ids or [])
    if operation:
        s["operation"] = F(operation, operation, pages)
    if temperature is not None:
        if isinstance(temperature, (int, float)):
            s["temperature"] = F(f"{temperature} °C", quantity(temperature, "°C", f"{temperature} °C"), pages)
        else:
            s["temperature"] = F(temperature, temperature, pages)
    if time is not None:
        if isinstance(time, tuple):
            number, unit = time
            s["time"] = F(f"{number} {unit}", quantity(number, unit, f"{number} {unit}"), pages)
        else:
            s["time"] = F(time, time, pages)
    if atmosphere:
        s["atmosphere"] = F(atmosphere, atmosphere, pages)
    if pressure:
        s["pressure"] = F(pressure, pressure, pages)
    r["steps"].append(s)


def add_outcome(r, oid, kind, targets, pages, *, raw, normalized=None, product_role=None,
                yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev(*pages))
    o["value"] = F(raw, normalized if normalized is not None else raw, pages)
    if product_role:
        o["product_role"] = F(product_role, product_role, pages)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, pages)
    if basis:
        o["measurement_basis"] = F(basis, basis, pages)
    if scope:
        o["scope"] = F(scope, scope, pages)
    r["outcomes"].append(o)


def set_procedure(r, text, pages):
    r["procedure"]["original_text"] = F(text, text, pages)


def set_classification(r, kinds, pages, limitation):
    """kinds is a sequence of (change_type, element1, element2, before, after)."""
    changes = []
    labels = {}
    for i, (kind, a, b, before, after) in enumerate(kinds, 1):
        change = {
            "change_id": f"b{i}",
            "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
            "before_order": before,
            "after_order": after,
            "change_type": kind,
        }
        changes.append(change)
        if kind in {"formed", "hydrogen_gain", "order_increased"}:
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(change["change_id"])
    derivation = "Partial net-change interpretation from the source reaction drawing; no complete atom mapping is claimed."
    r["classification"]["bond_changes"] = [
        st.fact(None, c, ev(*pages), method="structure_derivation", derivation=derivation) for c in changes
    ]
    r["classification"]["labels"] = st.fact(
        None, list(labels.values()), ev(*pages), method="structure_derivation", derivation=derivation)
    r["classification"]["coverage"] = st.fact(
        None, "partial", ev(*pages), method="inference", derivation=derivation)
    r["classification"]["limitations"] = [limitation, "Exact atom endpoints and complete hydrogen mapping require structure review."]


def class_annulation(r, pages, *, rings=2, hydrogen_shift=False):
    kinds = [("formed", "C", "C", 0, 1) for _ in range(rings)]
    if hydrogen_shift:
        kinds += [("hydrogen_loss", "C", "H", 1, 0), ("hydrogen_gain", "C", "H", 0, 1)]
    set_classification(r, kinds, pages, "Annulation/rearrangement classification is product-class level because exact atom mapping is unavailable.")


def class_none(r, pages, reason):
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = st.fact(reason, [], ev(*pages), method="inference", derivation=reason)
    r["classification"]["coverage"] = st.fact(reason, "not_applicable", ev(*pages), method="inference", derivation=reason)
    r["classification"]["limitations"] = [reason]


def class_unresolved(r, pages, reason):
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = st.fact(reason, [], ev(*pages), method="inference", derivation=reason)
    r["classification"]["coverage"] = st.fact(reason, "unresolved", ev(*pages), method="inference", derivation=reason)
    r["classification"]["limitations"] = [reason]


template = json.loads((RUN / "submitted.json").read_text())
molecules = []

# Generic prior-work structures and current-paper primary structures.
molecules += [
    generic_molecule(f"{DOC}:mol:generic_1", "1", (1,), "generic 1,7-allenediene skeleton shown in Scheme 1a", ["R", "R'", "X"], "Scheme 1a is a variable prior-work scaffold."),
    generic_molecule(f"{DOC}:mol:generic_2", "2", (1,), "generic 5,7-fused [4C+3C] annulation product", ["R", "R'", "X"], "Generic prior-work product; variables are not expanded."),
    generic_molecule(f"{DOC}:mol:generic_2prime", "2'", (1,), "generic stereoisomeric 5,7-fused [4C+3C] annulation product", ["R", "R'", "X"], "Distinct wedge/dash product drawing retained separately."),
    generic_molecule(f"{DOC}:mol:generic_3", "3", (1,), "generic 5,6-fused formal [4C+2C] annulation product", ["R", "R'", "X"], "Generic prior-work product; variables are not expanded."),
]

primary_specs = [
    ("1a", (1,2,3,4), "unresolved", "unspecified", None),
    ("1b", (2,3,4,5), "unresolved", "unspecified", None),
    ("Z-1b", (4,5), "unresolved", "specified", "cis-conjugated diene explicitly denoted Z-1b"),
    ("1c", (4,5), "unresolved", "unspecified", None),
    ("1d", (4,5), "unresolved", "unspecified", None),
    ("2a", (1,2), "specified", "unspecified", "relative wedge/dash stereochemistry drawn; product reported as not observed"),
    ("2b", (2,3,5), "specified", "unspecified", "relative wedge/dash stereochemistry drawn; product reported as not observed"),
    ("2d", (4,5), "unspecified", "unspecified", None),
    ("3a", (1,2,3,4), "specified", "unspecified", "relative wedge/dash stereochemistry depicted"),
    ("3b", (2,3,5), "specified", "unspecified", "relative wedge/dash stereochemistry depicted"),
    ("4a", (1,2,3,4), "specified", "unspecified", "relative stereochemistry drawn and stated to be confirmed by X-ray analysis"),
    ("4b", (2,3,4), "specified", "unspecified", "relative wedge/dash stereochemistry depicted"),
    ("4d", (4,5), "specified", "unspecified", "relative wedge/dash stereochemistry drawn; product not observed"),
    ("5a", (1,3,4), "specified", "unspecified", "relative wedge/dash stereochemistry depicted"),
    ("5b", (3,4,5,6), "specified", "unspecified", "relative stereochemistry depicted; desilylated derivative 5b' confirmed 6,6-bicyclic structure by X-ray"),
    ("5d", (4,5), "specified", "unspecified", "relative wedge/dash stereochemistry drawn; product not observed"),
    ("5b'", (3,), "specified", "unspecified", "desilylated X-ray derivative"),
    ("6b", (4,), "specified", "unspecified", "trans-fused product; complete diastereoselectivity stated"),
    ("6b'", (4,), "specified", "unspecified", "desilylated trans-fused product; complete diastereoselectivity stated"),
    ("7b", (4,), "specified", "unspecified", "relative wedge/dash stereochemistry depicted"),
    ("8b", (4,), "specified", "unspecified", "relative wedge/dash stereochemistry depicted"),
    ("9b", (4,5), "partially_specified", "unspecified", "cycloheptenone trapping product; stereochemical detail is incomplete in the main-text drawing"),
    ("10b", (4,), "partially_specified", "unspecified", "one alcohol product in the 10b/10b' pair"),
    ("10b'", (4,), "partially_specified", "unspecified", "second alcohol product in the 10b/10b' pair"),
]
for label, pages, chi, geo, assignment in primary_specs:
    mid = f"{DOC}:mol:{label.replace("'", '_prime').replace('-', '_').lower()}"
    molecules.append(molecule(mid, label, pages, chirality=chi, geometry=geo,
                              assignment=assignment))

def MID(label):
    if label in {"1", "2", "2'", "3"}:
        key = {"1":"generic_1","2":"generic_2","2'":"generic_2prime","3":"generic_3"}[label]
        return f"{DOC}:mol:{key}"
    return f"{DOC}:mol:{label.replace("'", '_prime').replace('-', '_').lower()}"


agents = {
    "Au1": "IPrAuCl/AgSbF6 (Au1)", "Au2'": "JohnPhosAuNTf2 (Au2')",
    "Au2": "JohnPhosAuCl/AgSbF6 (Au2)",
    "Au3": "(ArO)3PAuCl/AgSbF6 (Au3)", "Au4": "XPhosAuNTf2 (Au4)",
    "Au5": "RuPhosAuNTf2 (Au5)", "Au6": "IPrMeAuCl/AgNTf2 (Au6)",
    "Au7": "(C6F5)3PAuCl/AgNTf2 (Au7)", "AuCl3": "AuCl3", "AuBr3": "AuBr3",
    "PicAuCl2": "PicAuCl2", "PtCl2": "PtCl2",
    "PtEthene": "[PtCl2(ethene)2]2", "PtP": "PtCl2/P(C6F5)3P (1:1)",
    "CO": "carbon monoxide", "CH2Cl2": "CH2Cl2", "DCE": "1,2-DCE",
    "toluene": "toluene", "cyclohexane": "cyclohexane", "dioxane": "1,4-dioxane",
    "MeOH": "MeOH", "EtOH": "EtOH", "THF": "THF", "H2": "H2",
    "PdC": "Pd/C", "PdOH2C": "Pd(OH)2/C", "PtO2": "PtO2", "TBAF": "TBAF",
    "DMP": "Dess–Martin periodinane", "IPrAuNTf2": "IPrAuNTf2",
    "DCPO": "3,5-dichloropyridine oxide", "NaBH4": "NaBH4",
}
agent_pages = {
    "Au1": (1, 2, 3, 4), "Au2'": (2, 3), "Au2": (1,),
    "Au3": (1, 2, 3, 4), "Au4": (2, 3), "Au5": (2, 3),
    "Au6": (2, 3), "Au7": (1, 2, 3, 4), "AuCl3": (2, 3),
    "AuBr3": (2, 3), "PicAuCl2": (2, 3), "PtCl2": (1, 2, 3, 4, 5, 6),
    "PtEthene": (2, 3, 6), "PtP": (2, 3), "CO": (1, 3, 4, 5, 6),
    "CH2Cl2": (2, 3, 4), "DCE": (3,), "toluene": (2, 3),
    "cyclohexane": (3, 4, 5, 6), "dioxane": (3, 4),
    "MeOH": (4,), "EtOH": (4,), "THF": (4,), "H2": (4,),
    "PdC": (4,), "PdOH2C": (4,), "PtO2": (4,), "TBAF": (4,),
    "DMP": (4,), "IPrAuNTf2": (4, 5), "DCPO": (4,), "NaBH4": (4,),
}
for key, name in agents.items():
    pages = agent_pages[key]
    molecules.append(molecule(f"{DOC}:mol:agent_{key.lower()}", name, pages, name=name,
                              chirality="unresolved", geometry="unresolved",
                              note="Reagent/catalyst/solvent identity retained as printed; no external identity or stereochemical assertion added."))

def A(key): return f"{DOC}:mol:agent_{key.lower()}"

reactions = []

# Scheme 1a: two classes of cited prior work.
for suffix, product_labels, title, trans in [
    ("scheme1a_prior_4plus3", ["2", "2'"], "Previously developed divergent allenediene cycloadditions: formal [4C + 3C]", "formal [4C + 3C] annulation"),
    ("scheme1a_prior_4plus2", ["3"], "Previously developed divergent allenediene cycloadditions: formal [4C + 2C]", "formal [4C + 2C] annulation"),
]:
    r = reaction(suffix, (1,), record_type="literature_example", categories=["other"], name=title,
                 transformation=trans, location="Scheme 1a (cited prior work)")
    r["participants"].append(participant("substrate", MID("1"), "reactant", (1,)))
    catalyst_keys = ("PtCl2", "Au1", "Au2") if suffix.endswith("4plus3") else ("Au3",)
    for i, key in enumerate(catalyst_keys, 1):
        r["participants"].append(participant(f"catalyst_{i}", A(key), "catalyst", (1,)))
    for i, pl in enumerate(product_labels):
        r["participants"].append(participant(f"product_{i+1}", MID(pl), "product", (1,), product_role="prior-work product"))
    class_annulation(r, (1,), rings=2, hydrogen_shift=(suffix.endswith("4plus3")))
    set_procedure(r, "Scheme 1a summarizes previously reported catalyst-dependent allenediene cycloadditions.", (1,))
    reactions.append(r)

# Scheme 1b is a current-paper overview, not another set of specific experiments.
r = reaction("scheme1b_current_overview", (1,), record_type="generic_scheme", categories=["other"],
             name="Current divergent annulations", transformation="catalyst-dependent divergent annulation",
             location="Scheme 1b")
r["participants"] = [participant("substrate", MID("1a"), "reactant", (1,))]
for i, (key, role) in enumerate((("Au3","catalyst"),("Au7","catalyst"),("PtCl2","catalyst"),("Au1","catalyst"),("CO","reagent")), 1):
    r["participants"].append(participant(f"overview_agent_{i}", A(key), role, (1,)))
for pl in ("3a", "4a", "5a", "2a"):
    r["participants"].append(participant(f"product_{pl}", MID(pl), "product", (1,),
                                         product_role="not observed" if pl == "2a" else "overview product"))
class_annulation(r, (1,), rings=2, hydrogen_shift=True)
r["classification"]["limitations"].append("Overview contains divergent products with different detailed net changes; labels are a union-level summary.")
set_procedure(r, "Scheme 1b overview; exact conditions are represented in later specific records.", (1,))
reactions.append(r)

# Scheme 2 discovery experiment.
r = reaction("scheme2_1a_au1", (2,), categories=["other"], name="Unexpected reaction outcome of 1a using Au1",
             transformation="Au(I)-catalyzed allenediene annulation", location="Scheme 2")
r["participants"] = [
    participant("substrate", MID("1a"), "reactant", (2,)),
    participant("catalyst", A("Au1"), "catalyst", (2,), mol_percent=10, reference=MID("1a")),
    participant("p3a", MID("3a"), "product", (2,), product_role="minor product"),
    participant("p4a", MID("4a"), "product", (2,), product_role="isolated major product"),
    participant("p2a", MID("2a"), "product", (2,), product_role="not observed"),
]
set_solvent(r, A("CH2Cl2"), "CH2Cl2", (2,))
add_step(r, "reaction", (2,), participant_ids=["substrate", "catalyst"], operation="treated with Au1",
         temperature="rt", time=(2, "h"))
add_outcome(r, "ratio", "selectivity", [MID("3a"), MID("4a")], (2,), raw="3a : 4a = 1 : 2",
            normalized={"metric":"product ratio","order":["3a","4a"],"ratio":[1,2]}, scope="3a:4a")
add_outcome(r, "yield_4a", "yield", [MID("4a")], (2,), raw="57% yield",
            normalized=quantity(57, "%", "57%"), product_role="isolated major product",
            yield_type="isolated", basis="isolated yield")
add_outcome(r, "not_observed_2a", "no_reaction", [MID("2a")], (2,), raw="2a not observed", normalized="not observed")
class_annulation(r, (2,), rings=3, hydrogen_shift=True)
set_procedure(r, "1a; Au1 (10 mol%); CH2Cl2; rt; 2 h.", (2,))
reactions.append(r)

# Table 1, every row retained separately.
table1 = [
    (1,"Au1",10,"CH2Cl2","rt",2,"5 : 1",70,"4b",False),
    (2,"Au2'",10,"CH2Cl2","rt",12,None,None,None,True),
    (3,"Au4",10,"DCE",80,12,"1.1 : 1",18,"4b",False),
    (4,"Au5",10,"CH2Cl2",40,12,"16 : 1",32,"4b",False),
    (5,"Au6",10,"CH2Cl2","rt",2,"1 : 0",96,"4b",False),
    (6,"Au3",10,"CH2Cl2",-15,1.5,"1 : 7",66,"3b",False),
    (7,"Au7",5,"CH2Cl2",-15,1,"1 : 12",70,"3b",False),
    (8,"AuCl3",10,"CH2Cl2","rt",12,None,None,None,True),
    (9,"AuBr3",10,"CH2Cl2","rt",12,None,None,None,True),
    (10,"PicAuCl2",10,"toluene",110,8,"1 : 0",17,"4b",False),
    (11,"PtCl2",10,"toluene",110,0.5,"1 : 0",43,"4b",False),
    (12,"PtEthene",5,"toluene",110,1,"1 : 0",53,"4b",False),
    (13,"PtP",10,"toluene",110,1,"1 : 0",62,"4b",False),
    (14,"PtCl2",10,"toluene",70,2,"1 : 0",81,"4b",False),
    (15,"PtCl2",10,"cyclohexane",50,2,"1 : 0",88,"4b",False),
    (16,"PtCl2",10,"dioxane",30,2.5,"1 : 0",99,"4b",False),
    (17,"PtCl2",1,"dioxane",30,24,"1 : 0",92,"4b",False),
]
for entry, catalyst, loading, solvent, temp, hours, ratio, yld, yprod, complex_mix in table1:
    pages = (2,3) if entry <= 13 else (3,)
    r = reaction(f"table1_entry_{entry:02d}", pages, categories=["optimization"],
                 name="Exploration of the cycloaddition of 1b",
                 transformation="Au/Pt-catalyzed allenediene annulation", location=f"Table 1 entry {entry}")
    r["general_scheme_ids"] = [f"{DOC}:scheme1b_current_overview"]
    r["participants"] = [
        participant("substrate", MID("1b"), "reactant", pages),
        participant("catalyst", A(catalyst), "catalyst", pages, mol_percent=loading, reference=MID("1b")),
    ]
    # Rows marked only as a complex mixture do not identify 4b or 3b.  Keep
    # those molecules out of the participant list instead of turning table
    # header possibilities into observed products.  The bracketed 2b drawing
    # is explicitly marked "not observed" for the table and remains targeted.
    if not complex_mix:
        r["participants"].extend([
            participant("p4b", MID("4b"), "product", pages, product_role="ratio component"),
            participant("p3b", MID("3b"), "product", pages, product_role="ratio component"),
        ])
    r["participants"].append(
        participant("p2b", MID("2b"), "product", pages, product_role="not observed")
    )
    set_solvent(r, A(solvent), agents[solvent], pages, concentration=0.05, basis=MID("1b"))
    add_step(r, "reaction", pages, participant_ids=["substrate", "catalyst"], operation="cycloaddition screening",
             temperature=temp, time=(hours, "h"))
    if ratio:
        a, b = [x.strip() for x in ratio.split(":")]
        add_outcome(r, "ratio_4b_3b", "selectivity", [MID("4b"),MID("3b")], pages,
                    raw=f"4b : 3b = {ratio}", normalized={"metric":"4b:3b ratio","raw_ratio":ratio},
                    basis="1H NMR of the crude mixture", scope="4b:3b")
    if yld is not None:
        add_outcome(r, f"isolated_yield_{yprod}", "yield", [MID(yprod)], pages,
                    raw=f"{yld}% ({yprod})", normalized=quantity(yld, "%", f"{yld}%"),
                    product_role="isolated product", yield_type="isolated", basis="isolated yield")
    if complex_mix:
        add_outcome(r, "complex_mixture", "other", [], pages, raw="complex mixture of products",
                    normalized="complex mixture of products", scope="unidentified products")
        class_unresolved(r, pages, "Complex product mixture prevents product-specific net-change classification.")
    else:
        mixed = ratio is not None and not ratio.endswith(": 0")
        class_annulation(r, pages, rings=3 if yprod == "4b" else 2,
                         hydrogen_shift=(yprod == "4b" or mixed))
    add_outcome(r, "not_observed_2b", "no_reaction", [MID("2b")], pages,
                raw="2b not observed", normalized="not observed")
    set_procedure(r, f"Table 1 entry {entry}; conditions and footnotes preserved.", pages)
    reactions.append(r)

# Scheme 3: three CO-condition experiments.
scheme3 = [
    ("toluene",50,1,"8 : 1",59,False),
    ("cyclohexane",50,5,"8 : 1",79,False),
    ("cyclohexane",50,2,"9 : 1",85,True),
]
for i,(solvent,temp,hours,ratio,yld,wet) in enumerate(scheme3,1):
    r = reaction(f"scheme3_co_condition_{i}", (3,), categories=["optimization"],
                 name="PtCl2-catalyzed reaction of 1b under a CO atmosphere",
                 transformation="CO-dependent Pt-catalyzed annulation", location=f"Scheme 3 condition row {i}")
    r["participants"] = [
        participant("substrate", MID("1b"), "reactant", (3,)),
        participant("catalyst", A("PtCl2"), "catalyst", (3,), mol_percent=10, reference=MID("1b")),
        participant("co", A("CO"), "reagent", (3,)),
        participant("p5b", MID("5b"), "product", (3,), product_role="major product"),
        participant("p4b", MID("4b"), "product", (3,), product_role="minor product"),
    ]
    set_solvent(r,A(solvent),agents[solvent],(3,))
    add_step(r,"reaction",(3,),participant_ids=["substrate","catalyst","co"],operation="PtCl2-catalyzed annulation under CO",
             temperature=temp,time=(hours,"h"),atmosphere="CO (balloon)")
    add_outcome(r,"ratio","selectivity",[MID("5b"),MID("4b")],(3,),raw=f"5b:4b = {ratio}",
                normalized={"metric":"5b:4b ratio","raw_ratio":ratio},scope="5b:4b")
    add_outcome(r,"yield_5b","yield",[MID("5b")],(3,),raw=f"{yld}%",normalized=quantity(yld,"%",f"{yld}%"),product_role="major product")
    if wet:
        r["notes"].append(F("cyclohexane with 1 ppm content of water","cyclohexane with 1 ppm content of water",(3,)))
    class_annulation(r,(3,),rings=3,hydrogen_shift=True)
    set_procedure(r,f"PtCl2 (10 mol%)/CO; {agents[solvent]}; {temp} °C; {hours} h.",(3,))
    reactions.append(r)

# Scheme 4: three divergent outcomes from 1a; prose on page 4 repeats the figure.
scheme4 = [
    ("pt_4a","PtCl2",10,"cyclohexane",50,2,None,"4a",88,True),
    ("ptco_5a","PtCl2",10,"cyclohexane",50,2,"CO","5a",85,True),
    ("au7_3a","Au7",10,"CH2Cl2","rt",None,None,"3a",54,False),
]
for suffix,cat,load,solv,temp,hours,atm,prod,yld,hshift in scheme4:
    r = reaction(f"scheme4_{suffix}",(3,4),categories=["substrate_scope"],
                 name="Catalyst-dependent skeletal diversity from 1a",
                 transformation="divergent allenediene annulation",location=f"Scheme 4 arrow to {prod}")
    r["general_scheme_ids"]=[f"{DOC}:scheme1b_current_overview"]
    r["participants"]=[participant("substrate",MID("1a"),"reactant",(3,4)),
                       participant("catalyst",A(cat),"catalyst",(3,4),mol_percent=load,reference=MID("1a")),
                       participant("product",MID(prod),"product",(3,4),product_role="reported product")]
    if atm:
        r["participants"].append(participant("co",A("CO"),"reagent",(3,4)))
    set_solvent(r,A(solv),agents[solv],(3,4))
    add_step(r,"reaction",(3,4),participant_ids=[p["participant_id"] for p in r["participants"] if p["role"]!="product"],
             operation="catalyst-dependent annulation",temperature=temp,time=(hours,"h") if hours else None,
             atmosphere="CO" if atm else None)
    add_outcome(r,f"yield_{prod}","yield",[MID(prod)],(3,4),raw=f"{yld}% yield",normalized=quantity(yld,"%",f"{yld}%"),product_role="reported product")
    class_annulation(r,(3,4),rings=3 if prod in {"4a","5a"} else 2,hydrogen_shift=hshift)
    set_procedure(r,f"Scheme 4 conditions for {prod}.",(3,4))
    reactions.append(r)

# Scheme 5 derivatizations.
derivs = [
    ("scheme5_eq1a_6b","6b",[("H2","reagent",None),("PdC","catalyst",10),("MeOH","solvent",None)],"MeOH","rt",None,"H2 (1 atm)",99,"alkene hydrogenation with cyclopropane opening"),
    ("scheme5_eq1b_6bprime","6b'",[("H2","reagent",None),("PdOH2C","catalyst",10),("MeOH","solvent",None)],"MeOH","rt",None,"H2 (1 atm)",70,"hydrogenation, cyclopropane opening and TBS removal"),
    ("scheme5_eq2_7b","7b",[("H2","reagent",None),("PtO2","catalyst",10),("EtOH","solvent",None)],"EtOH","rt",None,"H2 (5 bar)",99,"selective double-bond hydrogenation"),
]
for suffix,prod,extra,solv,temp,hours,atm,yld,trans in derivs:
    r=reaction(suffix,(4,),categories=["derivatization"],name="Preliminary exploration of the synthetic potential of adduct 4b",transformation=trans,location="Scheme 5")
    r["participants"]=[participant("substrate",MID("4b"),"reactant",(4,)),participant("product",MID(prod),"product",(4,),product_role="reported product")]
    for key,role,load in extra:
        if role=="solvent": continue
        r["participants"].append(participant(key.lower(),A(key),role,(4,),mol_percent=load,reference=MID("4b")) if load else participant(key.lower(),A(key),role,(4,)))
    set_solvent(r,A(solv),agents[solv],(4,))
    pressure="5 bar" if "5 bar" in atm else "1 atm"
    add_step(r,"reaction",(4,),participant_ids=[p["participant_id"] for p in r["participants"] if p["role"]!="product"],operation=trans,temperature=temp,time=hours,atmosphere="H2",pressure=pressure)
    add_outcome(r,"yield","yield",[MID(prod)],(4,),raw=f"{yld}%",normalized=quantity(yld,"%",f"{yld}%"),product_role="reported product")
    if prod in {"6b","6b'"}:
        add_outcome(r,"diastereoselectivity","selectivity",[MID(prod)],(4,),raw="complete diastereoselectivity",normalized="complete diastereoselectivity",scope="diastereoselectivity")
    set_classification(r,[("order_decreased","C","C",2,1),("hydrogen_gain","C","H",0,1),("hydrogen_gain","C","H",0,1)]+([("broken","C","C",1,0)] if prod in {"6b","6b'"} else []),(4,),"Partial hydrogenation net changes; exact atoms are not mapped.")
    set_procedure(r,trans+"; duration not reported in main text.",(4,))
    reactions.append(r)

r=reaction("scheme5_eq3_8b",(4,),categories=["derivatization"],name="Conversion of the OTBS group of 4b to ketone 8b",transformation="desilylation and oxidation",location="Scheme 5 eq. 3")
r["participants"]=[participant("substrate",MID("4b"),"reactant",(4,)),participant("tbaf",A("TBAF"),"reagent",(4,)),participant("dmp",A("DMP"),"reagent",(4,)),participant("product",MID("8b"),"product",(4,),product_role="reported product")]
set_solvent(r,A("THF"),"THF",(4,))
r["participants"].append(participant("solvent_dcm",A("CH2Cl2"),"solvent",(4,)))
r["solvent_conditions"]["composition"]=F("THF then CH2Cl2",["THF","CH2Cl2"],(4,))
add_step(r,"desilylation",(4,),participant_ids=["substrate","tbaf"],operation="TBAF treatment",temperature="rt")
add_step(r,"oxidation",(4,),participant_ids=["substrate","dmp"],operation="Dess–Martin oxidation",temperature="rt")
add_outcome(r,"yield","yield",[MID("8b")],(4,),raw="79%",normalized=quantity(79,"%","79%"),product_role="reported ketone")
set_classification(r,[("broken","O","Si",1,0),("order_increased","C","O",1,2),("hydrogen_loss","C","H",1,0)],(4,),"Net desilylation/oxidation changes inferred from the source structures; atom mapping is incomplete.")
set_procedure(r,"1) TBAF, THF, rt; 2) Dess–Martin, CH2Cl2, rt; times and doses not reported.",(4,))
reactions.append(r)

# Scheme 6 mechanistic probes.
r=reaction("scheme6_eq1_4b_ptcl2co_control",(4,),categories=["control","mechanistic_study"],name="Test of 4b as a precursor to 5b",transformation="no reaction control",location="Scheme 6a eq. 1")
r["participants"]=[participant("substrate",MID("4b"),"reactant",(4,)),participant("catalyst",A("PtCl2"),"catalyst",(4,),mol_percent=10,reference=MID("4b")),participant("co",A("CO"),"reagent",(4,)),participant("expected_5b",MID("5b"),"product",(4,),product_role="not observed")]
set_solvent(r,A("cyclohexane"),"cyclohexane",(4,))
add_step(r,"reaction",(4,),participant_ids=["substrate","catalyst","co"],operation="standard PtCl2/CO conditions",temperature=50,atmosphere="CO (1 atm)")
add_outcome(r,"no_reaction","no_reaction",[MID("5b")],(4,),raw="no reaction; 5b not observed",normalized="no reaction")
class_none(r,(4,),"No identifiable transformation was observed; no expected-product bond class is assigned.")
set_procedure(r,"4b; PtCl2 (10 mol%); CO (1 atm); cyclohexane; 50 °C; time not reported in the supplied main text.",(4,))
reactions.append(r)

r=reaction("scheme6_eq2_oxidative_trapping",(4,5),categories=["mechanistic_study"],name="Oxidative trapping of the putative carbene",transformation="oxidative trapping followed by NaBH4 reduction",location="Scheme 6a eq. 2")
r["participants"]=[participant("substrate",MID("1b"),"reactant",(4,5)),participant("catalyst",A("IPrAuNTf2"),"catalyst",(4,5),mol_percent=10,reference=MID("1b")),participant("oxidant",A("DCPO"),"reagent",(4,5),equivalents=10,reference=MID("1b")),participant("nabh4",A("NaBH4"),"reagent",(4,5)),participant("p9b",MID("9b"),"product",(4,5),product_role="trapped ketone intermediate"),participant("p4b",MID("4b"),"product",(4,5),product_role="co-product"),participant("p10b",MID("10b"),"product",(4,5),product_role="reduction product"),participant("p10bp",MID("10b'"),"product",(4,5),product_role="reduction product")]
set_solvent(r,A("CH2Cl2"),"CH2Cl2 then MeOH",(4,5))
r["participants"].append(participant("solvent_meoh",A("MeOH"),"solvent",(4,5)))
r["solvent_conditions"]["composition"]=F("CH2Cl2 then MeOH",["CH2Cl2","MeOH"],(4,5))
add_step(r,"oxidative_trapping",(4,5),participant_ids=["substrate","catalyst","oxidant"],operation="IPrAuNTf2/3,5-dichloropyridine oxide treatment")
add_step(r,"reduction",(4,5),participant_ids=["p9b","nabh4"],operation="crude residue reduced with NaBH4 in MeOH")
add_outcome(r,"yield_4b","yield",[MID("4b")],(4,5),raw="67% yield",normalized=quantity(67,"%","67%"),product_role="co-product")
add_outcome(r,"approx_yield_9b","yield",[MID("9b")],(4,5),raw="approximately 30% yield",normalized=quantity(30,"%","~30%",qualifier="approx"),product_role="trapped ketone intermediate")
add_outcome(r,"overall_yield_10","yield",[MID("10b"),MID("10b'")],(4,5),raw="10b/10b' (30% overall yield)",normalized=quantity(30,"%","30%"),product_role="alcohol pair",scope="overall from 1b through trapping and reduction")
set_classification(r,[("formed","C","C",0,1),("formed","C","C",0,1),("formed","C","O",0,1),("hydrogen_gain","C","H",0,1),("formed","O","H",0,1)],(4,5),"Multi-stage union of oxidative trapping and reduction; product-specific atom mapping is unresolved.")
set_procedure(r,"1) IPrAuNTf2 (10 mol%), 3,5-dichloropyridine oxide (10 equiv), CH2Cl2; 2) NaBH4, MeOH. Temperatures/times not reported here.",(4,5))
reactions.append(r)

for label in ("Z-1b","1c"):
    r=reaction(f"scheme6_eq3_{label.lower().replace('-', '_')}_screen",(4,5),categories=["control","mechanistic_study"],name="Unreactive allenediene control",transformation="catalyst-screen control",location=f"Scheme 6a eq. 3 ({label})",unresolved=["The main text summarizes several catalysts in one arrow; catalyst-specific rows are only referenced to unavailable Tables S8–S9.","At 50 °C no product was reported; at higher temperature complex mixtures were reported; individual times and solvents are unavailable."])
    r["participants"]=[participant("substrate",MID(label),"reactant",(4,5))]
    for i,key in enumerate(("PtCl2","Au1","Au3"),1):
        r["participants"].append(participant(f"catalyst_{i}",A(key),"catalyst",(4,5),mol_percent=10,reference=MID(label)))
    add_step(r,"screen_at_50c",(4,5),participant_ids=[p["participant_id"] for p in r["participants"]],operation="screen with PtCl2, PtCl2/CO, Au1, Au3 and other catalysts",temperature=50)
    add_step(r,"higher_temperature_screen",(4,5),participant_ids=[p["participant_id"] for p in r["participants"]],operation="higher-temperature catalyst screen")
    add_outcome(r,"no_product_50c","no_reaction",[],(4,5),raw="no reaction @ 50 °C",normalized="no reaction at 50 °C")
    add_outcome(r,"complex_mixtures_high_temp","other",[],(4,5),raw="complex mixtures at high temperatures",normalized="complex mixtures at high temperatures")
    class_none(r,(4,5),"No defined product was obtained in the summarized controls; expected bond classes are not assigned.")
    set_procedure(r,"Scheme 6a eq. 3 catalyst-screen summary; exact row conditions require the unavailable SI.",(4,5))
    reactions.append(r)

r=reaction("scheme6b_preliminary_mechanistic_hypothesis",(4,5,6),record_type="proposed_mechanism",categories=["mechanistic_study"],name="Preliminary mechanistic hypothesis",transformation="proposed divergent carbene rearrangements",location="Scheme 6b")
r["participants"]=[participant("substrate",MID("1b"),"reactant",(4,5,6)),participant("p4b",MID("4b"),"product",(4,5,6),product_role="proposed C7-H insertion path"),participant("p3b",MID("3b"),"product",(4,5,6),product_role="proposed C10 alkyl migration path"),participant("p5b",MID("5b"),"product",(4,5,6),product_role="proposed C7 alkyl migration path")]
class_annulation(r,(4,5,6),rings=3,hydrogen_shift=True)
r["classification"]["limitations"].append("Arrow-pushing proposal is author interpretation, not an independently observed experiment.")
set_procedure(r,"Author-proposed [4C+3C] cycloaddition to carbene A' followed by divergent rearrangements.",(4,5,6))
reactions.append(r)

# Scheme 7, four distinct rows.
scheme7=[
    (1,"PtCl2","dioxane",30,12,None,99,False),
    (2,"PtCl2","cyclohexane",50,2,"CO",99,True),
    (3,"PtCl2","cyclohexane",50,3.5,None,83,True),
    (4,"Au1","CH2Cl2","rt",None,None,91,False),
]
for row,cat,solv,temp,hours,atm,yld,wet in scheme7:
    r=reaction(f"scheme7_row_{row}",(4,5),categories=["substrate_scope"],name="Pt(II)- and Au(I)-catalyzed annulations of 1d",transformation="formal [4C + 3C] annulation",location=f"Scheme 7 row {row}")
    r["participants"]=[participant("substrate",MID("1d"),"reactant",(4,5)),participant("catalyst",A(cat),"catalyst",(4,5),mol_percent=10,reference=MID("1d")),participant("p2d",MID("2d"),"product",(4,5),product_role="reported product"),participant("p4d",MID("4d"),"product",(4,5),product_role="not observed"),participant("p5d",MID("5d"),"product",(4,5),product_role="not observed")]
    if atm: r["participants"].append(participant("co",A("CO"),"reagent",(4,5)))
    set_solvent(r,A(solv),agents[solv],(4,5))
    add_step(r,"reaction",(4,5),participant_ids=[p["participant_id"] for p in r["participants"] if p["role"]!="product"],operation="annulation of 1d",temperature=temp,time=(hours,"h") if hours is not None else None,atmosphere="CO" if atm else None)
    add_outcome(r,"yield_2d","yield",[MID("2d")],(4,5),raw=f"{yld}%",normalized=quantity(yld,"%",f"{yld}%"),product_role="reported product")
    add_outcome(r,"not_observed_4d","no_reaction",[MID("4d")],(4,5),raw="4d not observed",normalized="not observed")
    add_outcome(r,"not_observed_5d","no_reaction",[MID("5d")],(4,5),raw="5d not observed",normalized="not observed")
    if wet: r["notes"].append(F("cyclohexane with 1 ppm content of water","cyclohexane with 1 ppm content of water",(4,5)))
    class_annulation(r,(4,5),rings=2,hydrogen_shift=True)
    set_procedure(r,f"Scheme 7 row {row}; catalyst 10 mol%; solvent/temperature/time as printed.",(4,5))
    reactions.append(r)

template["molecules"] = molecules
template["reactions"] = reactions
template["source_context"] = [
    {"context_id":f"{DOC}:context:fig1","kind":"structure_assignment","evidence_ids":ev(2),"description":"Fig. 1 reports an X-ray structure of 4a; natural-product structures are contextual, not synthesis records."},
    {"context_id":f"{DOC}:context:esi_notice","kind":"supplement_notice","evidence_ids":ev(1),"description":"The paper states that ESI exists, but no SI file was explicitly supplied or associated in this workspace."},
    {"context_id":f"{DOC}:context:dft_fig2_3","kind":"computational_analysis","evidence_ids":ev(5),"description":"DFT Gibbs-energy profiles in Figs. 2–3 are mechanistic calculations, not synthetic experiments."},
    {"context_id":f"{DOC}:context:dft_fig4","kind":"computational_analysis","evidence_ids":ev(6),"description":"Fig. 4 NBO orbitals, distances and Wiberg-index discussion are computational source context."},
]

st.write(RUN / "submitted.json", template)

inventory = {
    "document_id": DOC,
    "inspected_pages": [1,2,3,4,5,6],
    "visual_method": "full-page rendered JPEGs inspected at original detail",
    "candidate_structure_note": "A plausible connectivity-only 1b SMILES was graph-checked separately but not promoted because the exact tether connectivity and source stereochemistry were not sufficiently secure.",
    "regions": {
        "1": ["Scheme 1a prior-work [4C+3C] and [4C+2C] arrows", "Scheme 1b current-paper overview"],
        "2": ["Scheme 2 specific Au1 experiment", "Fig. 1 X-ray/context", "Table 1 prose lead-in"],
        "3": ["Table 1 entries 1–17 and footnotes", "Scheme 3 condition rows 1–3", "Scheme 4 arrows 1–3"],
        "4": ["Scheme 5 arrows 1a/1b/2/3", "Scheme 6 probe equations 1–3 and hypothesis", "Scheme 7 rows 1–4"],
        "5": ["Scheme 6/7 prose continuation", "DFT Figs. 2–3"],
        "6": ["DFT Fig. 4 and mechanistic discussion"],
    },
    "reaction_count": len(reactions),
    "molecule_count": len(molecules),
}
st.write(RUN / "raw_model/page_inventory.json", inventory)

ids = {r["reaction_id"] for r in reactions}
def rid(s):
    full=f"{DOC}:{s}"
    assert full in ids
    return full

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page":1,"state":"attempted","regions":[
            {"region":"Scheme 1a","finding":"Two cited prior-work annulation pathways were separated as literature examples.","evidence_ids":ev(1)},
            {"region":"Scheme 1b","finding":"Current-paper catalyst-divergence overview retained as one generic scheme; later arrows carry the specific conditions.","evidence_ids":ev(1)},
            {"region":"title/footnote","finding":"Title/DOI and ESI-availability statement read; ESI was not supplied and was not associated.","evidence_ids":ev(1)}],
         "record_ids":[rid("scheme1a_prior_4plus3"),rid("scheme1a_prior_4plus2"),rid("scheme1b_current_overview")],
         "unresolved":["Generic R/R'/X structures are not expanded to concrete molecules.","No supplied SI supports the experimental details referenced in the footnote."]},
        {"page":2,"state":"attempted","regions":[
            {"region":"Scheme 2","finding":"Specific 1a/Au1 experiment recorded with 3a:4a ratio, 4a isolated yield and 2a not observed.","evidence_ids":ev(2)},
            {"region":"Fig. 1","finding":"4a X-ray assignment and natural-product context retained as source context, not new syntheses.","evidence_ids":ev(2)},
            {"region":"right-column Table 1 discussion","finding":"Narrative descriptions of entries 1–13 were treated as duplicate support for rows enumerated on page 3.","evidence_ids":ev(2)}],
         "record_ids":[rid("scheme2_1a_au1")]+[rid(f"table1_entry_{i:02d}") for i in range(1,14)],
         "unresolved":["Exact isomeric SMILES for the fused products were not promoted from the drawing."]},
        {"page":3,"state":"attempted","regions":[
            {"region":"Table 1","finding":"All 17 optimization rows and all footnotes were enumerated separately, including three complex-mixture rows.","evidence_ids":ev(3)},
            {"region":"Scheme 3","finding":"Three PtCl2/CO condition rows were separated; the wet-cyclohexane footnote is attached only to row 3.","evidence_ids":ev(3)},
            {"region":"Scheme 4","finding":"Three catalyst-dependent arrows from 1a were separated as specific experiments.","evidence_ids":ev(3)}],
         "record_ids":[rid(f"table1_entry_{i:02d}") for i in range(1,18)]+[rid(f"scheme3_co_condition_{i}") for i in range(1,4)]+[rid("scheme4_pt_4a"),rid("scheme4_ptco_5a"),rid("scheme4_au7_3a")],
         "unresolved":["Yield basis for Schemes 3–4 is not stated in the supplied main text.","Exact structures of Au complexes are referred to unavailable Tables S1–S2."]},
        {"page":4,"state":"attempted","regions":[
            {"region":"Scheme 5","finding":"Four derivatization arrows were separated; the two eq. 1 catalyst variants are independent records.","evidence_ids":ev(4)},
            {"region":"Scheme 6a","finding":"No-reaction control, oxidative-trapping sequence, and two summarized substrate controls were separated.","evidence_ids":ev(4)},
            {"region":"Scheme 6b","finding":"Arrow-pushing hypothesis recorded as proposed mechanism, not a proved experiment.","evidence_ids":ev(4)},
            {"region":"Scheme 7","finding":"Four catalyst/solvent condition rows were separated; 4d and 5d remain not-observed outcomes.","evidence_ids":ev(4)}],
         "record_ids":[rid("scheme4_pt_4a"),rid("scheme4_ptco_5a"),rid("scheme4_au7_3a"),rid("scheme5_eq1a_6b"),rid("scheme5_eq1b_6bprime"),rid("scheme5_eq2_7b"),rid("scheme5_eq3_8b"),rid("scheme6_eq1_4b_ptcl2co_control"),rid("scheme6_eq2_oxidative_trapping"),rid("scheme6_eq3_z_1b_screen"),rid("scheme6_eq3_1c_screen"),rid("scheme6b_preliminary_mechanistic_hypothesis")]+[rid(f"scheme7_row_{i}") for i in range(1,5)],
         "unresolved":["Most Scheme 5–7 yield measurement bases and several durations/amounts are not reported.","Catalyst-specific conditions behind Scheme 6 eq. 3 require unavailable SI Tables S8–S9."]},
        {"page":5,"state":"attempted","regions":[
            {"region":"top/right prose","finding":"Continuation confirms Scheme 6 controls and Scheme 7 outcomes; linked to existing records rather than duplicated.","evidence_ids":ev(5)},
            {"region":"Figs. 2–3","finding":"DFT energy profiles and stationary-point structures retained as computational source context, not synthesis records.","evidence_ids":ev(5)}],
         "record_ids":[rid("scheme6_eq2_oxidative_trapping"),rid("scheme6_eq3_z_1b_screen"),rid("scheme6_eq3_1c_screen"),rid("scheme6b_preliminary_mechanistic_hypothesis")]+[rid(f"scheme7_row_{i}") for i in range(1,5)],
         "unresolved":["DFT structures were not converted into independent molecule/reaction experiments."]},
        {"page":6,"state":"attempted","regions":[
            {"region":"Fig. 4 and both text columns","finding":"NBO orbitals, distances, Wiberg indices and pathway discussion retained as computational/mechanistic context; no new experimental reaction arrow was present.","evidence_ids":ev(6)}],
         "record_ids":[rid("scheme6b_preliminary_mechanistic_hypothesis")],
         "unresolved":["Computational atom numbering was not used to assert complete experimental atom mappings."]},
    ],
}
st.write(RUN / "coverage/pages_0001_0006.json", coverage)

tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    page = int(task["evidence_id"].split(":page", 1)[1].split(":", 1)[0])
    if 1 <= page <= 6:
        task["state"] = "attempted"
        task["reason"] = "Original rendered page inspected; findings and record mappings are in coverage/pages_0001_0006.json."
st.write(RUN / "tasks.json", tasks)
print(json.dumps({"molecules":len(molecules),"reactions":len(reactions)},ensure_ascii=False))
