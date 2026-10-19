#!/usr/bin/env python3
"""Build the page-1--6 source-derived v3.1 submission for Wang et al., JACS 2021.

This file records the manual visual transcription used to create the machine
submission.  It deliberately leaves unreported procedures and unsafe dense
structures unresolved.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_252796c15b38495a"
EV = {n: f"{DOC}:page{n}:image" for n in range(1, 7)}
TXT = {n: f"{DOC}:page{n}:text" for n in range(1, 8)}


def q(number, unit, **extra):
    return {"number": number, "unit": unit, **extra}


def original(raw, normalized=None, ev=(), method="original_extraction"):
    return st.fact(raw, raw if normalized is None else normalized, ev, method=method)


def derived(normalized, ev, why):
    return st.fact(None, normalized, ev, method="structure_derivation", derivation=why)


def inferred(normalized, ev, why):
    return st.fact(None, normalized, ev, method="inference", derivation=why)


bundle = st.read(RUN / "submitted.json")
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
molecules = {}
stereo_diagnostics = {}


def add_molecule(mid, label, ev, *, name=None, smiles=None,
                 chirality="unresolved", geometry="unresolved",
                 representation=None, note=None, label_method="original_extraction"):
    if mid in molecules:
        return mid
    m = st.new_molecule(mid, DOC, list(ev))
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = st.fact(
        label if label_method == "original_extraction" else None,
        label,
        ev,
        method=label_method,
        derivation=("Compound label inferred from the product numbering relationship in Table 2."
                    if label_method == "inference" else None),
    )
    m["names"]["original"] = original(name or label, name or label, ev)
    m["names"]["english"] = original(name or label, name or label, ev)
    m["cas_numbers"] = st.fact()
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            smiles, ev,
            "Candidate connectivity transcribed from the displayed structure; RDKit parsing is a graph check, not source identity verification.",
        )
        diag = st.stereo_inspect(smiles)
        stereo_diagnostics[mid] = diag
        m["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], ev,
            "RDKit canonical isomeric SMILES generated from the source-derived candidate graph.",
        )
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], ev,
            "RDKit stereo graph key generated from the source-derived candidate graph.",
        )
        rep = representation or ("partial" if any(p["specified"] == "Unspecified" for p in diag["potential_stereo"]) else "not_applicable")
    else:
        m["representation_type"] = "unresolved"
        rep = representation or "not_encoded"
    m["stereochemistry"]["chirality_status"] = inferred(
        chirality, ev,
        "Status assigned from the displayed bonds/stereochemical marks; absence of a mark is retained as unspecified or unresolved, never racemic.",
    )
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, ev,
        "Status assigned from the displayed alkene/allene geometry and any explicit E/Z labels.",
    )
    m["stereochemistry"]["representation_status"] = inferred(
        rep, ev,
        "Compares source-visible stereochemistry with the candidate representation; unresolved drawings are not claimed to be encoded.",
    )
    m["identity_note"] = original(
        note or "Source label and drawing retained; no CAS number was reported and no external identity lookup was performed.",
        note or "Source label and drawing retained; no CAS number was reported and no external identity lookup was performed.",
        ev,
    )
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


def add_template(mid, label, ev, fixed, variables, note):
    m = st.new_molecule(mid, DOC, list(ev))
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = original(label, label, ev)
    m["names"]["original"] = original(label, label, ev)
    m["names"]["english"] = original(label, label, ev)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": original(fixed, fixed, ev),
        "variables": [{
            "label": v,
            "attachment_points": original("shown attachment point", "shown attachment point", ev),
            "definition": original(desc, desc, ev),
            "allowed_values": st.fact(),
            "constraints": st.fact(),
        } for v, desc in variables],
        "instance_mapping": [],
        "reason": original(note, note, ev),
    }
    m["stereochemistry"]["chirality_status"] = inferred("unresolved", ev, "Generic variables prevent a unique chirality assignment.")
    m["stereochemistry"]["geometry_status"] = inferred("unresolved", ev, "Generic variables prevent a unique geometry assignment.")
    m["stereochemistry"]["representation_status"] = inferred("not_encoded", ev, "The source is a generic R-group drawing, not a concrete stereochemical structure.")
    m["identity_note"] = original(note, note, ev)
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


def set_charge(p, key, raw, normalized, ev, ref=None, method="original_extraction"):
    p["charges"][key] = original(raw, normalized, ev, method=method)
    if key in ("equivalents", "mol_percent"):
        p["equivalence_reference"] = original(
            f"relative to {ref}", ref, ev, method=method
        )


def participant(pid, mid, role, ev, *, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, list(ev), stage=stage)
    if product_role:
        p["product_role"] = original(product_role, product_role, ev)
    return p


def classification(kind, ev):
    c = st.new_reaction(f"{DOC}:temporary", DOC, list(ev))["classification"]
    c["model_votes"] = []
    c["mapped_reaction_smiles"] = st.fact()
    c["limitations"] = ["Partial source-derived net-change description; no complete atom mapping or three-model vote was fabricated."]
    specs = {
        "core": [
            ("b1", "formed", "C", "C", 0, 1),
            ("b2", "broken", "C", "O", 1, 0),
            ("b3", "hydrogen_loss", "C", "H", 1, 0),
        ],
        "cyclization": [("b1", "formed", "C", "N", 0, 1)],
        "hydroboration_oxidation": [
            ("b1", "formed", "C", "O", 0, 1),
            ("b2", "formed", "O", "H", 0, 1),
            ("b3", "hydrogen_gain", "C", "H", 0, 1),
            ("b4", "order_decreased", "C", "C", 2, 1),
        ],
        "hydrosilylation": [
            ("b1", "formed", "C", "Si", 0, 1),
            ("b2", "hydrogen_gain", "C", "H", 0, 1),
            ("b3", "order_decreased", "C", "C", 2, 1),
        ],
        "substitution": [
            ("b1", "formed", "C", "S", 0, 1),
            ("b2", "broken", "C", "Cl", 1, 0),
        ],
        "complexation": [
            ("b1", "formed", "C", "Fe", 0, 1),
            ("b2", "formed", "C", "Fe", 0, 1),
        ],
        "exchange": [
            ("b1", "hydrogen_loss", "C", "H", 1, 0),
            ("b2", "hydrogen_gain", "C", "H", 0, 1),
        ],
    }
    if kind in ("unresolved", "not_applicable"):
        c["coverage"] = inferred(kind, ev, "No definite net bond-formation classification can be assigned to this conceptual or unresolved observation.")
        c["labels"] = inferred([], ev, "No bond-formation labels asserted.")
        c["bond_changes"] = []
        return c
    changes = []
    labels = {}
    for cid, ctype, a, b, before, after in specs[kind]:
        change = {
            "change_id": cid,
            "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
            "before_order": before,
            "after_order": after,
            "change_type": ctype,
        }
        if ctype.startswith("hydrogen_"):
            change["hydrogen_count"] = 1
        changes.append(derived(change, ev, "Net bond change inferred from the displayed reactant/product structures."))
        if ctype in ("formed", "hydrogen_gain", "order_increased"):
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(cid)
    c["bond_changes"] = changes
    c["labels"] = derived(list(labels.values()), ev, "Taxonomy labels calculated from the explicitly listed partial bond changes.")
    c["coverage"] = inferred("partial", ev, "Only source-visible net changes are represented; omitted coproducts and proton transfers prevent complete mapping.")
    return c


def reaction(rid, ev, *, record_type="specific_experiment", categories=(), page=1,
             region="", name=None, transform=None, class_kind="unresolved"):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, list(ev))
    r["record_type"] = record_type
    r["experiment_categories"] = original(list(categories), list(categories), ev)
    r["location"] = original(
        f"PDF page {page}, {region}",
        {"page": page, "region": region}, ev,
    )
    r["reaction_name"] = original(name, name, ev) if name else st.fact()
    if transform:
        r["transformation_type"] = inferred(transform, ev, "Concise transformation class inferred from the displayed reaction.")
    r["classification"] = classification(class_kind, ev)
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None,
        "checked_fields": [], "unresolved": [],
    }
    r["completeness"] = "partial"
    return r


def add_step(r, sid, ev, participant_ids, *, operation=None, temp=None, time=None,
             stage="reaction", atmosphere=None, pressure=None):
    s = st.new_step(sid, list(ev), stage=stage)
    s["participant_ids"] = list(participant_ids)
    if operation is not None:
        s["operation"] = original(operation, operation, ev)
    if temp is not None:
        s["temperature"] = original(temp[0], temp[1], ev)
    if time is not None:
        s["time"] = original(time[0], time[1], ev)
    if atmosphere is not None:
        s["atmosphere"] = original(atmosphere, atmosphere, ev)
    if pressure is not None:
        s["pressure"] = original(pressure, pressure, ev)
    r["steps"].append(s)


def add_outcome(r, oid, kind, targets, ev, raw, normalized, *, yield_type=None,
                basis=None, scope=None):
    o = st.new_outcome(oid, kind, list(targets), list(ev))
    o["value"] = original(raw, normalized, ev)
    if yield_type:
        o["yield_type"] = original(yield_type, yield_type, ev)
    if basis:
        o["measurement_basis"] = original(basis, basis, ev)
    if scope:
        o["scope"] = original(scope[0], scope[1], ev)
    r["outcomes"].append(o)


def add_note(r, text, ev, *, method="original_extraction", derivation=None):
    r["notes"].append(st.fact(
        text if method == "original_extraction" else None, text, ev,
        method=method, derivation=derivation,
    ))


# Concrete recurring molecules and generic Scheme 1 templates.
add_molecule("m_1a", "1a", [EV[2], EV[3]], name="monosubstituted undecatriene allene 1a (candidate)",
             smiles="C=C=CCCCCCCCC", chirality="achiral", geometry="not_applicable",
             note="The source shows terminal allene 1a bearing (CH2)7Me; the connectivity candidate preserves a terminal CH2 allene and no assigned stereochemistry.")
add_molecule("m_2a", "2a", [EV[2], EV[3], EV[4]], name="N-Mts,N-methyl hemiaminal methyl ether 2a",
             smiles="COC(c1ccc(Cl)cc1)N(C)S(=O)(=O)c1c(C)cc(C)cc1C",
             chirality="unspecified", geometry="not_applicable", representation="partial",
             note="The hemiaminal carbon is potentially stereogenic, but the source assigns no configuration.")
add_molecule("m_3a", "3a", [EV[2], EV[3], EV[4]], name="allenylic sulfonamide 3a",
             smiles="C=C=C(CCCCCCCC)C(c1ccc(Cl)cc1)N(C)S(=O)(=O)c1c(C)cc(C)cc1C",
             chirality="unspecified", geometry="not_applicable", representation="partial",
             note="The amino-bearing carbon is potentially stereogenic and is drawn without an absolute/relative assignment; the terminal CH2 allene is not treated as assigned axial chirality.")
add_template("m_scheme1_allene", "generic unactivated allene 1", [EV[1]], "H2C=C=CH-R1", [("R1", "alkyl substituent")], "Generic Scheme 1C substrate; R1 is not uniquely specified.")
add_template("m_scheme1_acetal", "generic N,O-acetal 2", [EV[1]], "Ar-CH(OMe)-N(R2)-PG", [("Ar", "aryl group"), ("R2", "N-alkyl group"), ("PG", "N-protecting group")], "Generic Scheme 1C iminium equivalent.")
add_template("m_scheme1_product", "generic allenylic amine product", [EV[1]], "H2C=C=C(R1)-CH(Ar)-N(R2)-PG", [("R1", "allene substituent"), ("Ar", "aryl group"), ("R2", "N-alkyl group"), ("PG", "N-protecting group")], "Generic 1,1-disubstituted allene product; configuration is not specified.")

# Common agents. Dense organometallic complexes remain unresolved rather than unsafe SMILES.
common = [
    ("m_cat_cpstar", "[Cp*Fe(CO)2(thf)]+BF4−", [EV[2]], "cationic Cp* iron catalyst"),
    ("m_cat_et5", "[(Et5C5)Fe(CO)2(thf)]+BF4−", [EV[2]], "cationic Et5C5 iron catalyst"),
    ("m_cat_npr5", "[(nPr5C5)Fe(CO)2(thf)]+BF4−", [EV[2], EV[3]], "cationic nPr5C5 iron catalyst"),
    ("m_bf3_etherate", "BF3·Et2O", [EV[2], EV[3], EV[4]], "boron trifluoride diethyl etherate"),
    ("m_26_lutidine", "2,6-lutidine", [EV[2]], "2,6-lutidine"),
    ("m_tmph", "TMPH", [EV[2]], "2,2,6,6-tetramethylpiperidine"),
    ("m_collidine", "2,4,6-collidine", [EV[2]], "2,4,6-collidine"),
    ("m_4cl_lutidine", "4-Cl-2,6-lutidine", [EV[2]], "4-chloro-2,6-lutidine"),
    ("m_4br_lutidine", "4-Br-2,6-lutidine", [EV[2], EV[3], EV[4]], "4-bromo-2,6-lutidine"),
    ("m_4br_picoline", "4-Br-2-picoline", [EV[2]], "4-bromo-2-picoline"),
    ("m_toluene", "toluene", [EV[2], EV[4]], "toluene"),
    ("m_phcf3", "PhCF3", [EV[2], EV[3], EV[4]], "benzotrifluoride"),
]
for mid, label, ev, name in common:
    add_molecule(mid, label, ev, name=name, chirality="achiral", geometry="not_applicable", representation="not_encoded")

# Page 1 inventory: conceptual/literature panels, current generic scheme, and proposed pathway.
for rid, rtype, region, name, transform in [
    ("scheme1A_literature_concept", "literature_example", "Scheme 1A", "Challenge: selective functionalization of monosubstituted allenes", "literature context"),
    ("scheme1B_prior_platform", "literature_example", "Scheme 1B", "Prior α-C−H functionalization of alkynes/alkenes", "prior-work α-C−H functionalization"),
]:
    r = reaction(rid, [EV[1]], record_type=rtype, page=1, region=region, name=name, transform=transform, class_kind="unresolved")
    r["review"]["unresolved"].append("Conceptual or cited prior-work panel; no current-paper experimental quantities or yield are asserted.")
    bundle["reactions"].append(r)

r = reaction("scheme1C_generic", [EV[1]], record_type="generic_scheme", page=1, region="Scheme 1C upper reaction arrow",
             name="α-Functionalization of unactivated alkylallenes (this report)", transform="allenic C(sp2)−H aminoalkylation", class_kind="core")
r["participants"] = [
    participant("p_allene", "m_scheme1_allene", "reactant", [EV[1]]),
    participant("p_acetal", "m_scheme1_acetal", "reactant", [EV[1]]),
    participant("p_product", "m_scheme1_product", "product", [EV[1]], product_role="generic primary product"),
]
r["review"]["unresolved"].append("Generic scheme only; R groups, charges, temperature and time are not specified here.")
bundle["reactions"].append(r)
r = reaction("scheme1C_proposed_pathway", [EV[1], TXT[2]], record_type="proposed_mechanism", page=1,
             region="Scheme 1C lower A/B/C/propargyliron pathway", name="Proposed iron-allene equilibration/deprotonation/SE2′ pathway",
             transform="proposed organoiron mechanism", class_kind="unresolved")
r["participants"] = [participant("p_allene", "m_scheme1_allene", "reactant", [EV[1]])]
r["review"]["unresolved"].append("Author-proposed intermediates A–C and propargyliron species are not treated as experimentally isolated products.")
bundle["reactions"].append(r)

# Table 1: 15 independent optimization rows.
pg_rows = {
    "CO2Me": ("m_2_co2me", "m_3_co2me"), "Bz": ("m_2_bz", "m_3_bz"),
    "Ts": ("m_2_ts", "m_3_ts"), "Ns": ("m_2_ns", "m_3_ns"),
    "Mbs": ("m_2_mbs", "m_3_mbs"), "Mts": ("m_2a", "m_3a"),
}
for pg, (smid, pmid) in pg_rows.items():
    if smid not in molecules:
        add_molecule(smid, f"2a-{pg}", [EV[2]], name=f"N,O-acetal with {pg} protecting group", chirality="unspecified", geometry="not_applicable", representation="not_encoded")
    if pmid not in molecules:
        add_molecule(pmid, f"3a-{pg}", [EV[2]], name=f"allenylic amine product with {pg} protecting group", chirality="unspecified", geometry="not_applicable", representation="not_encoded")

table1 = [
    (1, "Me", "m_26_lutidine", "CO2Me", 0), (2, "Me", "m_26_lutidine", "Bz", 0),
    (3, "Me", "m_26_lutidine", "Ts", 27), (4, "Me", "m_26_lutidine", "Ns", 21),
    (5, "Me", "m_26_lutidine", "Mbs", 32), (6, "Me", "m_26_lutidine", "Mts", 39),
    (7, "Me", "m_tmph", "Mts", 5), (8, "Me", "m_collidine", "Mts", 21),
    (9, "Me", "m_4cl_lutidine", "Mts", 50), (10, "Me", "m_4br_lutidine", "Mts", 55),
    (11, "Me", "m_4br_picoline", "Mts", 0), (12, "Me", "m_4br_lutidine", "Mts", 63),
    (13, "Et", "m_4br_lutidine", "Mts", 73), (14, "nPr", "m_4br_lutidine", "Mts", 77),
    (15, "nPr", "m_4br_lutidine", "Mts", 82),
]
cat_by_r = {"Me": "m_cat_cpstar", "Et": "m_cat_et5", "nPr": "m_cat_npr5"}
for entry, cat_r, base_mid, pg, yield_value in table1:
    smid, pmid = pg_rows[pg]
    r = reaction(f"table1_entry_{entry:02d}", [EV[2]], categories=["optimization"], page=2,
                 region=f"Table 1 entry {entry}", name="Optimization of allenylic amine synthesis",
                 transform="allenic C(sp2)−H aminoalkylation", class_kind="core")
    parts = [
        participant("p_1a", "m_1a", "reactant", [EV[2]]), participant("p_2", smid, "reactant", [EV[2]]),
        participant("p_cat", cat_by_r[cat_r], "catalyst", [EV[2]]), participant("p_bf3", "m_bf3_etherate", "reagent", [EV[2]]),
        participant("p_base", base_mid, "reagent", [EV[2]]), participant("p_solvent", "m_toluene" if entry <= 11 else "m_phcf3", "solvent", [EV[2]]),
        participant("p_product", pmid, "product", [EV[2]], product_role="target product"),
    ]
    r["participants"] = parts
    set_charge(parts[0], "equivalents", "2.0 equiv", q(2.0, "equiv"), [EV[2]], smid)
    set_charge(parts[1], "amount", "0.2 mmol scale (entry 15 isolated yield footnote)", q(0.2, "mmol"), [EV[2]]) if entry == 15 else None
    set_charge(parts[2], "mol_percent", "20 mol %", q(20, "mol%"), [EV[2]], smid)
    set_charge(parts[3], "equivalents", "2.0 equiv", q(2.0, "equiv"), [EV[2]], smid)
    base_eq = 4.0 if entry == 15 else 3.0
    set_charge(parts[4], "equivalents", f"{base_eq:.1f} equiv", q(base_eq, "equiv"), [EV[2]], smid)
    r["reaction_scale"] = original("2 (0.2 mmol) only explicitly printed for isolated entry 15", {"reference_molecule_id": smid, "amount": q(0.2, "mmol")}, [EV[2]]) if entry == 15 else st.fact()
    solvent = "toluene" if entry <= 11 else "PhCF3"
    r["solvent_conditions"]["composition"] = original(solvent, [solvent], [EV[2]])
    r["solvent_conditions"]["concentration"] = original("0.5 M", q(0.5, "M"), [EV[2]])
    r["solvent_conditions"]["concentration_basis"] = original("2a basis from displayed conditions", smid, [EV[2]])
    temp = 80 if entry <= 11 else 100
    add_step(r, "reaction", [EV[2]], [p["participant_id"] for p in parts[:-1]], operation="combine under Table 1 conditions",
             temp=(f"{temp} °C", q(temp, "°C")), time=("24 h", q(24, "h")))
    add_outcome(r, "yield_nmr", "yield", [pmid], [EV[2]], f"{yield_value}%", q(yield_value, "%"),
                yield_type="reported yield", basis="1H NMR using 1,1,2,2-tetrachloroethane as internal standard")
    if entry == 15:
        add_outcome(r, "yield_isolated", "yield", [pmid], [EV[2]], "79%", q(79, "%"),
                    yield_type="isolated yield", basis="isolated; 0.2 mmol scale")
    bundle["reactions"].append(r)

# Table 2 manual transcription.  Descriptions retain source-visible identity without unsafe full SMILES.
scope_top = [
    ("3a", 79, "n=6 linear alkyl", None), ("3b", 76, "n=8 linear alkyl", None),
    ("3c", 80, "cyclohexylethyl/β-branched alkyl", None), ("3d", 70, "cyclohexyl/secondary alkyl", None),
    ("3e", 65, "phenethyl-containing", None), ("3f", 71, "alkyl chloride-containing", None),
    ("3g", 69, "methyl ester-containing", None), ("3h", 64, "alkyl tosylate-containing", None),
    ("3i", 75, "4-bromophenoxy-containing", None), ("3j", 67, "2,6-dimethylphenoxy-containing", None),
    ("3k", 62, "aryl Bpin-containing", None), ("3l", 70, "nitrophenoxy-containing", None),
    ("3m", 40, "nitroaryl thioether-containing", "c"), ("3n", 67, "4-bromophenoxyethyl-containing", None),
    ("3o", 73, "acetate-protected alcohol", None), ("3p", 65, "benzoate-protected alcohol", None),
    ("3q", 71, "thiophene carboxylate-containing", None), ("3r", 68, "furan carboxylate-containing", None),
    ("3s", 73, "phthalimide-like ester-containing", None), ("3t", 61, "TBDMS ether-containing", None),
    ("3u", 73, "phthalimido-containing", None), ("3v", 55, "aryl alkyne-containing", None),
    ("3w", 77, "terminal alkene/dimethyl malonate-containing", None), ("3x", 34, "aryl allene-derived", None),
    ("3y", 69, "chloroxylenol-derived", None), ("3z", 62, "umbelliferone-derived", "c"),
    ("3aa", 64, "probenecid-derived", None), ("3ab", 78, "diflunisal-derived", "c"),
    ("3ac", 66, "gemfibrozil-derived", "c"), ("3ad", 71, "fenofibric-acid-derived", "c"),
]
scope_bottom = [
    ("4a", 65, "phenyl/N-methyl/Mts"), ("4b", 70, "4-fluorophenyl/N-methyl/Mts"),
    ("4c", 75, "4-bromophenyl/N-methyl/Mts"), ("4d", 77, "4-methoxycarbonylphenyl/N-methyl/Mts"),
    ("4e", 50, "4-methylphenyl/N-methyl/Mts"), ("4f", 75, "4-trifluoromethoxyphenyl/N-methyl/Mts"),
    ("4g", 69, "3-bromophenyl/N-methyl/Mts"), ("4h", 71, "benzodioxole-like difluoro/N-methyl/Mts"),
    ("4i", 52, "trifluoromethyl/methoxyphenyl/N-methyl/Mts"), ("4j", 61, "naphthyl/N-methyl/Mts"),
    ("4k", 88, "tert-butyl/N-methyl/Mts"), ("4l", 69, "4-chlorophenyl/N-ethyl/Mts"),
    ("4m", 58, "4-chlorophenyl/N-methyl/Ts"), ("4n", 66, "4-chlorophenyl/N-methyl/Mbs"),
    ("4o", 59, "phenyl/N-benzyl/Mbs"),
]


def scope_record(label, yield_value, desc, *, bottom=False, reverse=False, scale=None, cat15=False):
    ev = [EV[3], EV[4]]
    product_mid = "m_3a" if label == "3a" else f"m_{label}"
    add_molecule(product_mid, label, ev, name=f"allenylic sulfonamide {label}: {desc}",
                 chirality="unspecified", geometry="not_applicable", representation="not_encoded",
                 note="Concrete product drawing inspected; the amino-bearing carbon is shown without assigned configuration. Dense structure is retained by label/description rather than an unsafe complete SMILES.")
    if bottom:
        substrate_mid = "m_1a"
        reagent_mid = f"m_2_{label}"
        add_molecule(reagent_mid, f"2-series precursor to {label}", ev,
                     name=f"N,O-acetal precursor corresponding to {label}: {desc}", chirality="unspecified",
                     geometry="not_applicable", representation="not_encoded", label_method="inference",
                     note="Precursor identity is linked to the displayed product panel and N,O-acetal scope heading; exact precursor label is inferred and needs review.")
    else:
        substrate_mid = f"m_1_{label}"
        reagent_mid = "m_2a"
        add_molecule(substrate_mid, f"1-series precursor to {label}", ev,
                     name=f"allene precursor corresponding to {label}: {desc}", chirality="unspecified",
                     geometry="not_applicable", representation="not_encoded", label_method="inference",
                     note="Precursor identity is reconstructed only at the label/functional-group level from the displayed product panel; exact full SMILES is not asserted.")
    suffix = f"_{scale['tag']}" if scale else ""
    cats = ["substrate_scope"] + (["scale_up"] if scale else [])
    r = reaction(f"table2_{label}{suffix}", ev, categories=cats, page=3,
                 region=f"Table 2 product panel {label}" + (f" ({scale['raw']})" if scale else ""),
                 name="Scope of the allenylic amine synthesis", transform="allenic C(sp2)−H aminoalkylation", class_kind="core")
    parts = [
        participant("p_allene", substrate_mid, "reactant", ev), participant("p_acetal", reagent_mid, "reactant", ev),
        participant("p_cat", "m_cat_npr5", "catalyst", ev), participant("p_bf3", "m_bf3_etherate", "reagent", ev),
        participant("p_base", "m_4br_lutidine", "reagent", ev), participant("p_solvent", "m_phcf3", "solvent", ev),
        participant("p_product", product_mid, "product", ev, product_role="target product"),
    ]
    r["participants"] = parts
    ref_mid = substrate_mid if reverse else reagent_mid
    if scale:
        set_charge(parts[0] if reverse else parts[1], "amount", scale["raw"], q(scale["mmol"], "mmol"), [EV[3]])
        r["reaction_scale"] = original(scale["raw"], {"reference_molecule_id": ref_mid, "amount": q(scale["mmol"], "mmol")}, [EV[3]])
    elif reverse:
        set_charge(parts[0], "amount", "1 (0.2 mmol)", q(0.2, "mmol"), [EV[3]])
        r["reaction_scale"] = original("1 (0.2 mmol)", {"reference_molecule_id": substrate_mid, "amount": q(0.2, "mmol")}, [EV[3]])
    else:
        set_charge(parts[1], "amount", "2 (0.2 mmol)", q(0.2, "mmol"), [EV[3]])
        r["reaction_scale"] = original("2 (0.2 mmol)", {"reference_molecule_id": reagent_mid, "amount": q(0.2, "mmol")}, [EV[3]])
    if reverse:
        set_charge(parts[1], "equivalents", "2 (1.5 equiv)", q(1.5, "equiv"), [EV[3]], ref_mid)
    else:
        set_charge(parts[0], "equivalents", "1 (2.0 equiv)", q(2.0, "equiv"), [EV[3]], ref_mid)
    set_charge(parts[2], "mol_percent", f"{15 if cat15 else 20} mol %", q(15 if cat15 else 20, "mol%"), [EV[3]], ref_mid)
    set_charge(parts[3], "equivalents", "BF3·Et2O (2.0 equiv)", q(2.0, "equiv"), [EV[3]], ref_mid)
    set_charge(parts[4], "equivalents", "4-bromolutidine (4.0 equiv)", q(4.0, "equiv"), [EV[3]], ref_mid)
    r["solvent_conditions"]["composition"] = original("PhCF3", ["PhCF3"], [EV[3]])
    r["solvent_conditions"]["concentration"] = original("0.5 M", q(0.5, "M"), [EV[3]])
    r["solvent_conditions"]["concentration_basis"] = original("limiting substrate basis", ref_mid, [EV[3]])
    add_step(r, "reaction", [EV[3]], [p["participant_id"] for p in parts[:-1]], operation="combine under Table 2 conditions",
             temp=("100 °C", q(100, "°C")))
    add_outcome(r, "yield", "yield", [product_mid], [EV[3]], f"{yield_value}%", q(yield_value, "%"), yield_type="reported yield")
    add_outcome(r, "regioselectivity", "selectivity", [product_mid], [EV[4]], ">20:1 over other regioisomers",
                {"qualifier": "gt", "ratio": [20, 1], "order": [label, "other regioisomers"]},
                basis="1H NMR analysis of the crude material",
                scope=(">20:1 regioisomer ratio", {"metric": "regioisomer ratio", "order": [label, "other regioisomers"], "sample": "crude"}))
    if label in ("3v", "3w") and not scale:
        competitor = "alkyne" if label == "3v" else "alkene"
        add_outcome(r, "chemoselectivity", "selectivity", [product_mid], [EV[2]], f">20:1 allene over {competitor} functionalization",
                    {"qualifier": "gt", "ratio": [20, 1], "order": ["allene functionalization", f"{competitor} functionalization"]},
                    basis="1H NMR analysis of the crude material",
                    scope=(f">20:1 allene/{competitor}", {"metric": "chemoselectivity", "order": ["allene", competitor], "sample": "crude"}))
    r["review"]["unresolved"].extend([
        "Table 2 does not report reaction time or atmosphere.",
        "Yield measurement basis is not printed in the supplied main-text table and remains unreported.",
        "Dense structures are retained by source label/description unless a safe candidate graph was available.",
    ])
    bundle["reactions"].append(r)


for label, y, desc, foot in scope_top:
    scope_record(label, y, desc, reverse=(foot == "c"))
for label, y, desc in scope_bottom:
    scope_record(label, y, desc, bottom=True)
scope_record("3a", 71, "n=6 linear alkyl", scale={"tag": "scale3mmol", "raw": "3.0 mmol scale", "mmol": 3.0}, cat15=True)
scope_record("3ad", 65, "fenofibric-acid-derived", scale={"tag": "scale5mmol", "raw": "2.348 g, 5.0 mmol scale", "mmol": 5.0}, cat15=True)
scope_record("4k", 85, "tert-butyl/N-methyl/Mts", bottom=True, scale={"tag": "scale1p5mmol", "raw": "1.5 mmol scale", "mmol": 1.5}, cat15=True)

# Scheme 2 products/reagents.
scheme2_mols = [
    ("m_4a", "4a", "allenylic sulfonamide 4a", "unspecified", "not_applicable"),
    ("m_3f", "3f", "alkyl chloride-containing allenylic sulfonamide 3f", "unspecified", "not_applicable"),
    ("m_5", "5", "1,2,3-trisubstituted 3-pyrroline 5", "unspecified", "unspecified"),
    ("m_6z", "(Z)-6", "major Z allylic alcohol 6", "unspecified", "specified"),
    ("m_6e", "(E)-6", "minor E allylic alcohol 6", "unspecified", "specified"),
    ("m_7z", "(Z)-7", "major Z allylsilane 7", "unspecified", "specified"),
    ("m_7e", "(E)-7", "minor E allylsilane 7", "unspecified", "specified"),
    ("m_8", "8", "tetrazole-containing product 8", "unspecified", "not_applicable"),
]
for mid, lab, nm, ch, ge in scheme2_mols:
    if mid not in molecules:
        add_molecule(mid, lab, [EV[4]], name=nm, chirality=ch, geometry=ge, representation="not_encoded",
                     note="Scheme 2 structure inspected; any shown Z/E descriptor is retained textually because the dense full graph is not safely encoded.")
for mid, descriptor in [("m_6z", "Z"), ("m_6e", "E"), ("m_7z", "Z"), ("m_7e", "E")]:
    molecules[mid]["stereochemistry"]["double_bonds"] = original(
        f"({descriptor}) label in Scheme 2", [{"configuration": descriptor, "location": "new allylic alkene"}], [EV[4]]
    )
    molecules[mid]["stereochemistry"]["assignment"] = original(
        f"product labelled ({descriptor})-{mid.split('_')[1][0]}",
        {"configuration": descriptor, "assignment_scope": "new allylic alkene"}, [EV[4]]
    )
for mid, lab, name in [
    ("m_na_naph", "Na/C10H8", "sodium naphthalenide"), ("m_cuotf2", "Cu(OTf)2", "copper(II) triflate"),
    ("m_9bbn", "9-BBN", "9-borabicyclo[3.3.1]nonane"), ("m_h2o2", "H2O2", "hydrogen peroxide"),
    ("m_naoh", "NaOH", "sodium hydroxide"), ("m_tetrazole_thiol", "1-phenyl-1H-tetrazole-5-thiol", "1-phenyl-1H-tetrazole-5-thiol"),
    ("m_k2co3", "K2CO3", "potassium carbonate"), ("m_ki", "KI", "potassium iodide"),
    ("m_moco6", "Mo(CO)6", "molybdenum hexacarbonyl"), ("m_ph2sih2", "Ph2SiH2", "diphenylsilane"),
    ("m_thf", "THF", "tetrahydrofuran"), ("m_ch2cl2", "CH2Cl2", "dichloromethane"),
    ("m_dmf", "DMF", "N,N-dimethylformamide"),
]:
    add_molecule(mid, lab, [EV[4]], name=name, chirality="achiral", geometry="not_applicable", representation="not_encoded")

# 4a -> 5, two-step overall.
r = reaction("scheme2_product5", [EV[4]], categories=["derivatization"], page=4, region="Scheme 2, arrow a–b",
             name="Reductive desulfonylation followed by copper-catalyzed cyclization", transform="desulfonylation/cyclization", class_kind="cyclization")
r["participants"] = [participant("p_start", "m_4a", "reactant", [EV[4]]), participant("p_na", "m_na_naph", "reagent", [EV[4]]),
                     participant("p_thf", "m_thf", "solvent", [EV[4]]), participant("p_cu", "m_cuotf2", "catalyst", [EV[4]]),
                     participant("p_dcm", "m_ch2cl2", "solvent", [EV[4]]), participant("p_product", "m_5", "product", [EV[4]], product_role="target product")]
add_step(r, "desulfonylation", [EV[4]], ["p_start", "p_na", "p_thf"], operation="Na/C10H8 in THF", temp=("−78 °C", q(-78, "°C")), time=("8 h", q(8, "h")))
add_step(r, "cyclization", [EV[4]], ["p_start", "p_cu", "p_dcm"], operation="Cu(OTf)2 in CH2Cl2", temp=("rt", "rt"), time=("24 h", q(24, "h")))
add_outcome(r, "overall_yield", "yield", ["m_5"], [EV[4]], "51%", q(51, "%"), yield_type="two-step overall yield",
            scope=("4a to 5 across steps a–b", {"start": "4a", "end": "5", "steps": ["a", "b"]}))
r["review"]["unresolved"].append("Yield measurement basis is not stated in Scheme 2.")
bundle["reactions"].append(r)

# 3a -> 6, two-step hydroboration/oxidation.
r = reaction("scheme2_product6", [EV[4]], categories=["derivatization"], page=4, region="Scheme 2, arrow c–d",
             name="Hydroboration–oxidation", transform="allene hydroboration–oxidation to allylic alcohol", class_kind="hydroboration_oxidation")
r["participants"] = [participant("p_start", "m_3a", "reactant", [EV[4]]), participant("p_9bbn", "m_9bbn", "reagent", [EV[4]]),
                     participant("p_thf", "m_thf", "solvent", [EV[4]]), participant("p_h2o2", "m_h2o2", "reagent", [EV[4]]),
                     participant("p_naoh", "m_naoh", "reagent", [EV[4]]), participant("p_z", "m_6z", "product", [EV[4]], product_role="major product"),
                     participant("p_e", "m_6e", "product", [EV[4]], product_role="minor geometric isomer")]
set_charge(r["participants"][4], "concentration" if False else "volume", None, None, [EV[4]]) if False else None
add_step(r, "hydroboration", [EV[4]], ["p_start", "p_9bbn", "p_thf"], operation="9-BBN in THF", temp=("rt", "rt"), time=("5 h", q(5, "h")))
add_step(r, "oxidation", [EV[4]], ["p_h2o2", "p_naoh"], operation="H2O2, NaOH (3 M)", temp=("0 °C", q(0, "°C")), time=("2 h", q(2, "h")))
add_outcome(r, "yield", "yield", ["m_6z", "m_6e"], [EV[4]], "75%", q(75, "%"), yield_type="combined product yield")
add_outcome(r, "z_e", "selectivity", ["m_6z", "m_6e"], [EV[4]], "9:1 Z/E, crude", {"ratio": [9, 1], "order": ["Z", "E"]},
            scope=("9:1 Z/E, crude", {"metric": "Z/E", "order": ["m_6z", "m_6e"], "sample": "crude"}))
r["review"]["unresolved"].extend(["Yield measurement basis is not stated.", "The Z/E analytical method is not stated; only crude-sample scope is printed."])
bundle["reactions"].append(r)

# 3f -> 8 substitution.
r = reaction("scheme2_product8", [EV[4]], categories=["derivatization"], page=4, region="Scheme 2, arrow e",
             name="Nucleophilic substitution to tetrazole-containing product", transform="alkyl chloride substitution", class_kind="substitution")
r["participants"] = [participant("p_start", "m_3f", "reactant", [EV[4]]), participant("p_thiol", "m_tetrazole_thiol", "reactant", [EV[4]]),
                     participant("p_k2co3", "m_k2co3", "reagent", [EV[4]]), participant("p_ki", "m_ki", "additive", [EV[4]]),
                     participant("p_dmf", "m_dmf", "solvent", [EV[4]]), participant("p_product", "m_8", "product", [EV[4]], product_role="target product")]
add_step(r, "reaction", [EV[4]], [p["participant_id"] for p in r["participants"][:-1]], operation="substitution in DMF",
         temp=("90 °C", q(90, "°C")), time=("24 h", q(24, "h")))
add_outcome(r, "yield", "yield", ["m_8"], [EV[4]], "82%", q(82, "%"), yield_type="reported yield")
r["review"]["unresolved"].append("Yield measurement basis and reagent quantities are not stated in Scheme 2.")
bundle["reactions"].append(r)

# 3a -> 7 hydrosilylation.
r = reaction("scheme2_product7", [EV[4]], categories=["derivatization"], page=4, region="Scheme 2, arrow f",
             name="Molybdenum-catalyzed hydrosilylation", transform="allene hydrosilylation", class_kind="hydrosilylation")
r["participants"] = [participant("p_start", "m_3a", "reactant", [EV[4]]), participant("p_mo", "m_moco6", "catalyst", [EV[4]]),
                     participant("p_silane", "m_ph2sih2", "reagent", [EV[4]]), participant("p_solvent", "m_toluene", "solvent", [EV[4]]),
                     participant("p_z", "m_7z", "product", [EV[4]], product_role="major product"), participant("p_e", "m_7e", "product", [EV[4]], product_role="minor geometric isomer")]
add_step(r, "reaction", [EV[4]], [p["participant_id"] for p in r["participants"][:-2]], operation="Mo(CO)6, Ph2SiH2 in PhMe",
         temp=("110 °C", q(110, "°C")), time=("36 h", q(36, "h")))
add_outcome(r, "yield", "yield", ["m_7z", "m_7e"], [EV[4]], "80%", q(80, "%"), yield_type="combined product yield")
add_outcome(r, "z_e", "selectivity", ["m_7z", "m_7e"], [EV[4]], ">20:1 Z/E, crude", {"qualifier": "gt", "ratio": [20, 1], "order": ["Z", "E"]},
            scope=(">20:1 Z/E, crude", {"metric": "Z/E", "order": ["m_7z", "m_7e"], "sample": "crude"}))
r["review"]["unresolved"].extend(["Yield measurement basis is not stated.", "The Z/E analytical method is not stated; only crude-sample scope is printed."])
bundle["reactions"].append(r)

# Scheme 3 and prose mechanistic experiments.
mechanistic_smiles = {
    "m_1n": "C=C=CCCOc1ccc(Br)cc1",
    "m_1b": "C=C=CCCCCCCCCCC",
    "m_phccme": "CC#Cc1ccccc1",
}
for mid, lab, name, ch, ge in [
    ("m_1n", "1n", "4-bromophenoxyethyl-substituted terminal allene 1n", "achiral", "not_applicable"),
    ("m_fpstar_i", "Fp*I", "Fp*I", "unresolved", "not_applicable"),
    ("m_agbf4", "AgBF4", "silver tetrafluoroborate", "achiral", "not_applicable"),
    ("m_1n_complex", "[1n·Fp*]+BF4−", "cationic Fp*–allene complex", "specified", "not_applicable"),
    ("m_1n_d", "1n-d", "α-deuterated allene 1n-d", "achiral", "not_applicable"),
    ("m_nd4so4", "(ND4)2SO4", "ammonium-d4 sulfate", "achiral", "not_applicable"),
    ("m_nai", "NaI", "sodium iodide", "achiral", "not_applicable"),
    ("m_acetone", "acetone", "acetone", "achiral", "not_applicable"),
    ("m_alkyne_complex", "[η2-(PhC≡CMe)Fp*]+BF4−", "cationic Fp*–alkyne complex", "unresolved", "not_applicable"),
    ("m_1b", "1b", "long-chain terminal allene 1b", "achiral", "not_applicable"),
    ("m_1b_complex", "[1b·Fp*]+BF4−", "cationic Fp*–allene complex of 1b", "unresolved", "not_applicable"),
    ("m_phccme", "PhC≡CMe", "1-phenyl-1-propyne", "achiral", "not_applicable"),
    ("m_fpstar_thf", "[Fp*(thf)]+BF4−", "cationic Fp*–THF complex", "unresolved", "not_applicable"),
    ("m_cd2cl2", "CD2Cl2", "dichloromethane-d2", "achiral", "not_applicable"),
]:
    smi = mechanistic_smiles.get(mid)
    add_molecule(mid, lab, [EV[4]], name=name, smiles=smi, chirality=ch, geometry=ge,
                 representation=("partial" if smi and "C=C=C" in smi else "not_applicable" if smi else "not_encoded"),
                 note=("The source assigns the isolated allene complex as the anti regio-/stereoisomer by NMR and X-ray; organometallic stereochemistry is retained textually."
                       if mid == "m_1n_complex" else None))
molecules["m_1n_complex"]["stereochemistry"]["assignment"] = original(
    "single major anti regio- and stereoisomer (>20:1 r.r. and d.r.); anti confirmed by X-ray",
    {"relative_configuration": "anti", "rr": ">20:1", "dr": ">20:1", "assignment_method": ["1H NMR", "X-ray crystallography"]},
    [EV[4], TXT[4]],
)

r = reaction("scheme3A_complex_formation", [EV[4]], categories=["mechanistic_study"], page=4, region="Scheme 3A upper arrow",
             name="Stoichiometric preparation of [1n·Fp*]+BF4−", transform="allene–iron complex formation", class_kind="complexation")
r["participants"] = [participant("p_1n", "m_1n", "reactant", [EV[4]]), participant("p_fpi", "m_fpstar_i", "reactant", [EV[4]]),
                     participant("p_ag", "m_agbf4", "reagent", [EV[4]]), participant("p_tol", "m_toluene", "solvent", [EV[4]]),
                     participant("p_product", "m_1n_complex", "product", [EV[4]], product_role="isolated organoiron complex")]
add_step(r, "complexation", [EV[4]], ["p_1n", "p_fpi", "p_ag", "p_tol"], operation="combine 1n and Fp*I in the presence of AgBF4",
         temp=("50 °C", q(50, "°C")))
add_outcome(r, "yield", "yield", ["m_1n_complex"], [EV[4]], "78%", q(78, "%"), yield_type="reported yield")
add_outcome(r, "rr_dr", "selectivity", ["m_1n_complex"], [TXT[4]], ">20:1 r.r. and d.r.", {"qualifier": "gt", "ratio": [20, 1], "metrics": ["r.r.", "d.r."]},
            basis="1H NMR (rt to −55 °C)", scope=("single major anti regio- and stereoisomer", {"assignment": "anti", "metrics": ["r.r.", "d.r."]}))
r["review"]["unresolved"].extend(["Reaction time and yield measurement basis are not printed.", "Organometallic η2 bonding and anti stereochemistry are not encoded in SMILES."])
bundle["reactions"].append(r)

r = reaction("scheme3A_complex_reactivity", [EV[4]], categories=["mechanistic_study"], page=4, region="Scheme 3A lower arrow",
             name="Catalytic-conditions reactivity of [1n·Fp*]+BF4−", transform="allenic C(sp2)−H aminoalkylation from isolated complex", class_kind="core")
r["participants"] = [participant("p_complex", "m_1n_complex", "reactant", [EV[4]]), participant("p_acetal", "m_2a", "reactant", [EV[4]]),
                     participant("p_bf3", "m_bf3_etherate", "reagent", [EV[4]]), participant("p_base", "m_4br_lutidine", "reagent", [EV[4]]),
                     participant("p_solvent", "m_phcf3", "solvent", [EV[4]]), participant("p_product", "m_3n", "product", [EV[4]], product_role="target product")]
add_step(r, "reaction", [EV[4]], [p["participant_id"] for p in r["participants"][:-1]], operation="N,O-acetal (2a), BF3·Et2O, 4-bromolutidine in PhCF3",
         temp=("100 °C", q(100, "°C")))
add_outcome(r, "yield", "yield", ["m_3n"], [EV[4]], "39% NMR yield", q(39, "%"), yield_type="NMR yield", basis="NMR; internal standard not stated here")
r["review"]["unresolved"].append("Reagent quantities and reaction time are not printed in Scheme 3A.")
bundle["reactions"].append(r)

r = reaction("scheme3B_hd_exchange", [EV[4]], categories=["mechanistic_study"], page=4, region="Scheme 3B",
             name="H+/D+ exchange experiment", transform="allenic α-C(sp2) H/D exchange", class_kind="exchange")
r["participants"] = [participant("p_complex", "m_1n_complex", "reactant", [EV[4]]), participant("p_dsource", "m_nd4so4", "reagent", [EV[4]]),
                     participant("p_base", "m_4br_lutidine", "reagent", [EV[4]]), participant("p_solvent", "m_phcf3", "solvent", [EV[4]]),
                     participant("p_nai", "m_nai", "reagent", [EV[4]]), participant("p_acetone", "m_acetone", "solvent", [EV[4]]),
                     participant("p_1n", "m_1n", "product", [EV[4]], product_role="unlabeled component"), participant("p_1nd", "m_1n_d", "product", [EV[4]], product_role="deuterated component")]
set_charge(r["participants"][1], "equivalents", "(ND4)2SO4 (1.2 equiv)", q(1.2, "equiv"), [EV[4]], "m_1n_complex")
set_charge(r["participants"][2], "equivalents", "4-bromolutidine (1.5 equiv)", q(1.5, "equiv"), [EV[4]], "m_1n_complex")
add_step(r, "exchange", [EV[4]], ["p_complex", "p_dsource", "p_base", "p_solvent"], operation="expose complex to D+ source and base",
         temp=("80 °C", q(80, "°C")), time=("24 h", q(24, "h")))
add_step(r, "decomplexation", [EV[4]], ["p_nai", "p_acetone"], operation="then NaI in acetone", temp=("rt", "rt"))
add_outcome(r, "d_incorporation", "other", ["m_1n_d"], [EV[4]], "16% D", q(16, "%"),
            scope=("16% D exclusively at allenic α-C(sp2)", {"metric": "D incorporation", "site": "allenic α-C(sp2)"}))
r["review"]["unresolved"].append("D-incorporation analytical method and decomplexation time are not printed in Scheme 3B.")
bundle["reactions"].append(r)

r = reaction("scheme3C_alkyne_displacement", [EV[4]], categories=["mechanistic_study"], page=4, region="Scheme 3C",
             name="Displacement of coordinated alkyne by allene", transform="ligand displacement at iron", class_kind="complexation")
r["participants"] = [participant("p_alkyne_complex", "m_alkyne_complex", "reactant", [EV[4]]), participant("p_1b", "m_1b", "reactant", [EV[4]]),
                     participant("p_solvent", "m_cd2cl2", "solvent", [EV[4]]), participant("p_allene_complex", "m_1b_complex", "product", [EV[4]], product_role="allene complex"),
                     participant("p_alkyne", "m_phccme", "product", [EV[4]], product_role="displaced alkyne")]
set_charge(r["participants"][1], "equivalents", "1b (1.0 equiv)", q(1.0, "equiv"), [EV[4]], "m_alkyne_complex")
add_step(r, "displacement", [EV[4]], ["p_alkyne_complex", "p_1b", "p_solvent"], operation="ligand displacement in CD2Cl2",
         temp=("32 °C", q(32, "°C")), time=("15 h", q(15, "h")))
add_outcome(r, "conversion", "conversion", ["m_alkyne_complex"], [EV[4]], "91% conversion", q(91, "%"), basis="NMR context in scheme; exact conversion method label is not separately printed")
add_outcome(r, "allene_complex_yield", "yield", ["m_1b_complex"], [EV[4]], "85% NMR yield", q(85, "%"), yield_type="NMR yield", basis="NMR")
add_outcome(r, "alkyne_yield", "yield", ["m_phccme"], [EV[4]], "91% NMR yield", q(91, "%"), yield_type="NMR yield", basis="NMR")
r["review"]["unresolved"].append("The scheme does not name an internal standard or detailed NMR quantitation procedure.")
bundle["reactions"].append(r)

r = reaction("mechanism_competition_timecourse", [TXT[4]], categories=["mechanistic_study"], page=4, region="right-column prose below Scheme 3C",
             name="Competition between alkyne and allene complexation", transform="competitive ligand binding/time-course equilibration", class_kind="complexation")
r["participants"] = [participant("p_fp", "m_fpstar_thf", "reactant", [TXT[4]]), participant("p_1b", "m_1b", "reactant", [TXT[4]]),
                     participant("p_alkyne", "m_phccme", "reactant", [TXT[4]]), participant("p_solvent", "m_cd2cl2", "solvent", [TXT[4]]),
                     participant("p_allene_complex", "m_1b_complex", "product", [TXT[4]], product_role="allene complex"), participant("p_alkyne_complex", "m_alkyne_complex", "product", [TXT[4]], product_role="alkyne complex")]
set_charge(r["participants"][1], "equivalents", "allene 1b (2.0 equiv)", q(2.0, "equiv"), [TXT[4]], "m_fpstar_thf")
set_charge(r["participants"][2], "equivalents", "PhC≡CMe (2.0 equiv)", q(2.0, "equiv"), [TXT[4]], "m_fpstar_thf")
add_step(r, "early_timepoint", [TXT[4]], [p["participant_id"] for p in r["participants"][:4]], operation="competitive complexation in CD2Cl2",
         temp=("23 °C", q(23, "°C")), time=("10 min", q(10, "min")))
add_step(r, "equilibrated_timepoint", [TXT[4]], [p["participant_id"] for p in r["participants"][:4]], operation="prolonged equilibration",
         temp=("23 °C", q(23, "°C")), time=("68 h", q(68, "h")))
add_outcome(r, "early_ratio", "selectivity", ["m_1b_complex", "m_alkyne_complex"], [TXT[4]], "~1:2.5 allene/alkyne complexes at 10 min", {"qualifier": "approx", "ratio": [1, 2.5], "order": ["allene complex", "alkyne complex"], "time": "10 min"},
            scope=("10 min, <5% conversion", {"metric": "complex ratio", "time": q(10, "min"), "conversion": {"qualifier": "lt", "number": 5, "unit": "%"}}))
add_outcome(r, "early_conversion", "conversion", ["m_fpstar_thf"], [TXT[4]], "<5% conversion at 10 min", {"qualifier": "lt", "number": 5, "unit": "%", "raw_text": "<5%"})
add_outcome(r, "late_ratio", "selectivity", ["m_1b_complex", "m_alkyne_complex"], [TXT[4]], ">20:1 allene/alkyne complexes after 68 h", {"qualifier": "gt", "ratio": [20, 1], "order": ["allene complex", "alkyne complex"], "time": "68 h"},
            scope=("68 h, full conversion", {"metric": "complex ratio", "time": q(68, "h"), "conversion": "full"}))
add_outcome(r, "late_conversion", "conversion", ["m_fpstar_thf"], [TXT[4]], "full conversion after 68 h", "full conversion")
r["review"]["unresolved"].append("The prose does not state the analytical measurement basis for the complex ratios/conversions.")
bundle["reactions"].append(r)

r = reaction("mechanism_propargyl_deprotonation_attempt", [TXT[4]], categories=["mechanistic_study", "control"], page=4,
             region="left-column prose below Scheme 2", name="Attempted deprotonation/isolation of a propargyliron complex",
             transform="attempted organoiron deprotonation", class_kind="unresolved")
r["participants"] = [participant("p_complex", "m_1n_complex", "reactant", [TXT[4]])]
add_step(r, "attempt", [TXT[4]], ["p_complex"], operation="attempted deprotonation of [1n·Fp*]+BF4−",
         time=("over several hours", "several hours"))
add_outcome(r, "complex_mixture", "decomposition", [], [TXT[4]], "slow formation of a complex mixture of iron-containing species", "complex mixture")
r["review"]["unresolved"].extend(["The base, temperature and exact time for the current failed isolation attempt are not stated in the supplied prose.", "No definitive product identity or bond classification is possible."])
bundle["reactions"].append(r)

# Non-synthetic source context, including characterization and no-experiment pages.
bundle["source_context"] = [
    {"context_id": f"{DOC}:page1_metadata", "document_id": DOC, "evidence_ids": [EV[1], TXT[1]], "kind": "article_metadata",
     "value": {"title": "Iron-Catalyzed Contrasteric Functionalization of Allenic C(sp2)−H Bonds: Synthesis of α-Aminoalkyl 1,1-Disubstituted Allenes", "doi": "10.1021/jacs.1c07512", "journal": "J. Am. Chem. Soc. 2021, 143, 14998–15004"}},
    {"context_id": f"{DOC}:complex_characterization", "document_id": DOC, "evidence_ids": [EV[4], TXT[4]], "kind": "characterization",
     "value": "[1n·Fp*]+BF4− was a single major regio-/stereoisomer (>20:1 r.r. and d.r.) by 1H NMR from rt to −55 °C; X-ray assigned anti and an allene angle of 149.7°. This characterization is linked to the complex-formation experiment, not duplicated as a synthesis."},
    {"context_id": f"{DOC}:si_notice", "document_id": DOC, "evidence_ids": [EV[5], TXT[5]], "kind": "associated_content",
     "value": "The article states that Supporting Information contains experimental procedures and spectra, but no SI file was supplied or associated in this workspace."},
    {"context_id": f"{DOC}:page5_references", "document_id": DOC, "evidence_ids": [EV[5]], "kind": "no_current_experiment",
     "value": "Page 5 contains associated-content/author metadata and references; cited reactions are not current-paper experiments."},
    {"context_id": f"{DOC}:page6_references", "document_id": DOC, "evidence_ids": [EV[6]], "kind": "no_current_experiment",
     "value": "Page 6 contains references only; no current-paper reaction arrow, table row, or experimental observation was identified."},
]

# Raw visual transcription is separate from the formal bundle.
visual = {
    "document_id": DOC,
    "claimed_pages": [1, 2, 3, 4, 5, 6],
    "table1": [{"entry": e, "catalyst_R": cr, "base": b, "PG": pg, "yield_percent": y} for e, cr, b, pg, y in table1],
    "table2_top": [{"label": l, "yield_percent": y, "description": d, "footnote": f} for l, y, d, f in scope_top],
    "table2_bottom": [{"label": l, "yield_percent": y, "description": d} for l, y, d in scope_bottom],
    "scale_up": [{"label": "3a", "yield_percent": 71, "scale": "3.0 mmol", "footnote": "b"},
                 {"label": "3ad", "yield_percent": 65, "mass": "2.348 g", "scale": "5.0 mmol", "footnote": "b"},
                 {"label": "4k", "yield_percent": 85, "scale": "1.5 mmol", "footnote": "b"}],
    "scheme2": ["4a→5 (a–b), 51%", "3a→6 (c–d), 75%, 9:1 Z/E crude", "3f→8 (e), 82%", "3a→7 (f), 80%, >20:1 Z/E crude"],
    "scheme3": ["1n→[1n·Fp*]+BF4−, 78%", "isolated complex→3n, 39% NMR", "H/D exchange, 16% D", "alkyne displacement, 91% conversion/85% and 91% NMR yields"],
}

# Coverage records map every inspected page to reaction IDs or explicit no-experiment findings.
page_records = {p: [] for p in range(1, 7)}
for r in bundle["reactions"]:
    loc = r["location"]["normalized_value"]
    if isinstance(loc, dict) and loc.get("page") in page_records:
        page_records[loc["page"]].append(r["reaction_id"])
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/abstract", "finding": "Article metadata and abstract; no discrete experiment.", "evidence_ids": [EV[1], TXT[1]]},
            {"region": "Scheme 1A", "finding": "Conceptual literature challenge, classified as literature context.", "evidence_ids": [EV[1]]},
            {"region": "Scheme 1B", "finding": "Cited prior-work platform, not a current-paper experiment.", "evidence_ids": [EV[1]]},
            {"region": "Scheme 1C", "finding": "Current generic transformation plus author-proposed organoiron pathway; kept as generic scheme and proposed mechanism.", "evidence_ids": [EV[1]]},
        ], "record_ids": page_records[1], "unresolved": ["Generic R groups and proposed intermediates are not concrete experiments or fully specified structures."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Table 1 rows 1–15 and footnotes a–d", "finding": "Fifteen independent optimization experiments; NMR and isolated yields kept separate for entry 15.", "evidence_ids": [EV[2], TXT[2]]},
            {"region": "scope discussion", "finding": "Narrative cross-check for functional-group scope and 3v/3w chemoselectivity; no extra table rows invented.", "evidence_ids": [EV[2], TXT[2]]},
        ], "record_ids": page_records[2], "unresolved": ["Table 1 does not print the 0.2 mmol scale for entries 1–14; only entry 15 footnote explicitly states it."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Table 2 top: 3a–3ad", "finding": "Thirty displayed allene-scope products transcribed; three marked alternative-condition entries and two independent scale-up lines retained.", "evidence_ids": [EV[3]]},
            {"region": "Table 2 bottom: 4a–4o", "finding": "Fifteen N,O-acetal-scope products plus the independent 4k scale-up line transcribed.", "evidence_ids": [EV[3]]},
            {"region": "Table 2 footnotes a–c", "finding": "Default charges, 15 mol% catalyst exceptions, and reversed 1/2 stoichiometry encoded explicitly.", "evidence_ids": [EV[3], TXT[3]]},
        ], "record_ids": page_records[3], "unresolved": ["Reaction time, atmosphere and yield measurement basis are not printed.", "Many dense structures remain label/description-level and need expert structure review."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 2 arrows a–f", "finding": "Four independent derivatization experiments; two-step overall yields and Z/E outcomes preserved.", "evidence_ids": [EV[4], TXT[4]]},
            {"region": "Scheme 3A", "finding": "Stoichiometric complex preparation and a separate complex-reactivity experiment.", "evidence_ids": [EV[4], TXT[4]]},
            {"region": "Scheme 3B", "finding": "H/D exchange mechanistic experiment with 16% D incorporation.", "evidence_ids": [EV[4], TXT[4]]},
            {"region": "Scheme 3C and adjacent prose", "finding": "Alkyne displacement, competitive time-course, and unsuccessful propargyliron isolation attempt recorded separately.", "evidence_ids": [EV[4], TXT[4]]},
            {"region": "NMR/X-ray discussion", "finding": "Characterization kept in source_context and linked to complex formation, not duplicated as synthesis.", "evidence_ids": [EV[4], TXT[4]]},
        ], "record_ids": page_records[4], "unresolved": ["Several mechanistic arrows omit exact quantities, times, or analytical methods.", "Organometallic bonding/stereochemistry is not fully encoded."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "Associated Content", "finding": "SI availability is stated, but no SI was supplied; association remains unestablished.", "evidence_ids": [EV[5], TXT[5]]},
            {"region": "author information/acknowledgments/references", "finding": "No current-paper reaction experiment; cited work not imported as new results.", "evidence_ids": [EV[5]]},
        ], "record_ids": page_records[5], "unresolved": ["Experimental procedures and spectra reside in unavailable SI."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "references 4–11", "finding": "Reference list only; no current-paper table, reaction arrow, spectrum, or experiment.", "evidence_ids": [EV[6], TXT[6]]},
        ], "record_ids": page_records[6], "unresolved": []},
    ],
}

(RUN / "raw_model/visual_transcription.json").write_text(json.dumps(visual, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/stereo_diagnostics.json").write_text(json.dumps(stereo_diagnostics, ensure_ascii=False, indent=2) + "\n")
(RUN / "coverage/pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
(lambda tasks: (RUN / "tasks.json").write_text(json.dumps([
    ({**task, "state": "attempted", "reason": "Original page image inspected in the page-1–6 extraction session; see coverage/pages_0001_0006.json."}
     if task.get("document_id") == DOC and any(task.get("evidence_id") == TXT[p] for p in range(1, 7)) else task)
    for task in tasks
], ensure_ascii=False, indent=2) + "\n"))(json.loads((RUN / "tasks.json").read_text()))
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "page_counts": {str(k): len(v) for k, v in page_records.items()}}, ensure_ascii=False, indent=2))
