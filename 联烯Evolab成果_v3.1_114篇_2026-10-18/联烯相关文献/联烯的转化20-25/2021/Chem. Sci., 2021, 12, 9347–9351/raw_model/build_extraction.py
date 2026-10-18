#!/usr/bin/env python3
"""Build the source-derived v3.1 bundle for doc_3252f1d5f252557b.

This is a deterministic transcription of the visual page inventory.  It does
not certify source or chemical accuracy; merge quarantines every record as
needs_review.
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
DATASET = RUN / "dataset.json"
DOC = "doc_3252f1d5f252557b"
EV = {p: f"{DOC}:page{p}:image" for p in range(1, 6)}
EVT = {p: f"{DOC}:page{p}:text" for p in range(1, 6)}


def dump(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def F(raw=None, norm=None, ev=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, norm, ev, method=method, derivation=derivation,
                   status=status)


def inferred(norm, ev, why, raw=None):
    return F(raw, norm, ev, method="inference", derivation=why)


def derived(norm, ev, why, raw=None):
    return F(raw, norm, ev, method="structure_derivation", derivation=why)


data = st.read(DATASET)
assert not data.get("experiments"), "Refuse to build over existing experiments"
assert "standard" not in data, "Refuse to build over an existing v3 bundle"
doc = next(d for d in data["documents"] if d["document_id"] == DOC)
doc["title"] = "A palladium-catalyzed approach to allenic aromatic ethers and first total synthesis of terricollene A"
st.write(DATASET, data)

bundle = st.empty(data["documents"])
molecules = {}


def mol(mid, label, page, name=None, smiles=None, stereo="unresolved"):
    """Register one source molecule/agent with explicit stereo state."""
    full = f"{DOC}:mol_{mid}"
    if full in molecules:
        return full
    ev = [EV[page]]
    m = st.new_molecule(full, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    original = name or label
    m["names"]["original"] = F(original, original, ev)
    m["names"]["english"] = F(original, original, ev)
    m["names"]["abbreviations"] = F()
    m["names"]["aliases"] = F()
    m["cas_numbers"] = F()
    m["identity_note"] = F(
        "Identity retained from the source label/drawing; no CAS or external database enrichment was supplied.",
        "Source-labelled identity; structure remains a machine candidate unless explicitly encoded below.", ev)
    why = "Stereo state assessed from the supplied page drawing/name; absence of a descriptor is not treated as racemic."
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smiles, ev, "Candidate SMILES transcribed from the compact source drawing; parsing is not identity verification.", smiles)
        m["generic_template"] = None
        m["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], ev,
            "RDKit canonical isomeric SMILES calculated from the machine candidate; source accuracy remains unverified.")
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], ev,
            "Stereo graph key calculated from the machine candidate; it is not a sample identity.")
        m["stereochemistry"]["atom_centers"] = derived(diag["atom_centers"], ev, "RDKit graph diagnostic on the candidate SMILES.")
        m["stereochemistry"]["double_bonds"] = derived(diag["double_bonds"], ev, "RDKit graph diagnostic on the candidate SMILES.")
        m["stereochemistry"]["axial"] = F()
        m["stereochemistry"]["mixture"] = F()
        m["stereochemistry"]["assignment"] = F()
        if any(p["specified"] == "Unspecified" for p in diag["potential_stereo"]):
            # RDKit reports an unassigned potential element for the cumulene
            # candidates.  Retain that uncertainty instead of declaring the
            # source compound achiral or fully encoded.
            m["stereochemistry"]["chirality_status"] = inferred("unspecified", ev, "Toolkit diagnostic reports a potential unassigned cumulene/double-bond element; the source prints no assignment.")
            m["stereochemistry"]["geometry_status"] = inferred("unspecified", ev, "Toolkit diagnostic reports an unassigned potential element; no E/Z label is printed.")
            m["stereochemistry"]["representation_status"] = inferred("partial", ev, "Connectivity is encoded, but the toolkit flags an unspecified potential stereochemical element.")
        else:
            m["stereochemistry"]["chirality_status"] = inferred("achiral", ev, "Candidate graph and source drawing contain no potential stereogenic center or substituted allene axis.")
            m["stereochemistry"]["geometry_status"] = inferred("not_applicable", ev, "Candidate graph contains no independently configurable E/Z double bond.")
            m["stereochemistry"]["representation_status"] = inferred("encoded", ev, "All stereochemical information applicable to this candidate graph is encoded; source identity still needs review.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F()
        m["generic_template"] = None
        m["stereochemistry"]["atom_centers"] = F()
        m["stereochemistry"]["double_bonds"] = F()
        m["stereochemistry"]["axial"] = F()
        m["stereochemistry"]["mixture"] = F()
        m["stereochemistry"]["assignment"] = F()
        m["stereochemistry"]["chirality_status"] = inferred("unresolved", ev, why)
        m["stereochemistry"]["geometry_status"] = inferred("unresolved", ev, why)
        m["stereochemistry"]["representation_status"] = inferred("unresolved", ev, "Compact drawing was not converted into a complete stereochemically faithful SMILES.")
        m["stereochemistry"]["isomeric_smiles"] = F()
        m["stereochemistry"]["identity_key"] = F()
    if stereo == "axial_unspecified":
        m["stereochemistry"]["axial"] = inferred(
            [{"axis": "allene", "assignment": "unspecified"}], ev,
            "The 1,3-disubstituted allene drawing can support axial chirality, but no R_a/S_a or P/M assignment is printed.",
            "substituted allene drawn without an axial descriptor")
        m["stereochemistry"]["chirality_status"] = inferred("unspecified", ev, why)
        m["stereochemistry"]["geometry_status"] = inferred("unspecified", ev, "No separate alkene E/Z assignment is printed for the cumulene drawing.")
        m["stereochemistry"]["representation_status"] = inferred("not_encoded", ev, "Axial configuration is not assigned and no full stereochemical SMILES is asserted.")
    elif stereo == "central_specified":
        m["stereochemistry"]["atom_centers"] = F(
            "Wedge/dash stereochemistry is shown in the source drawing; individual R/S descriptors are not printed.",
            [{"assignment": "specified by drawing", "absolute_labels": "unreported"}], ev)
        m["stereochemistry"]["chirality_status"] = inferred("specified", ev, why)
        m["stereochemistry"]["geometry_status"] = inferred("not_applicable", ev, "No separate E/Z element is shown for this molecule.")
        m["stereochemistry"]["representation_status"] = inferred("not_encoded", ev, "Source stereocenters are retained textually; no stereochemical SMILES is asserted.")
    elif stereo == "geometry_Z":
        m["stereochemistry"]["double_bonds"] = F("(Z)-7", [{"assignment": "Z", "source_label": "(Z)-7"}], ev)
        m["stereochemistry"]["chirality_status"] = inferred("unresolved", ev, why)
        m["stereochemistry"]["geometry_status"] = inferred("specified", ev, why)
        m["stereochemistry"]["representation_status"] = inferred("not_encoded", ev, "Printed Z assignment retained, but full structure was not encoded.")
    elif stereo == "geometry_mixture":
        m["stereochemistry"]["double_bonds"] = F("E/Z mixture reported for this product", [{"assignment": "mixture", "ratio_in_reaction_outcome": True}], ev)
        m["stereochemistry"]["mixture"] = F("E/Z mixture", "E/Z mixture", ev)
        m["stereochemistry"]["chirality_status"] = inferred("unresolved", ev, why)
        m["stereochemistry"]["geometry_status"] = inferred("mixture", ev, why)
        m["stereochemistry"]["representation_status"] = inferred("not_encoded", ev, "Mixture identity and ratio are retained without inventing separate isomer structures.")
    elif stereo == "geometry_drawn":
        m["stereochemistry"]["double_bonds"] = F("Alkene geometry is depicted in the source drawing without an E/Z text label.", [{"assignment": "specified by drawing", "absolute_label": "unreported"}], ev)
        m["stereochemistry"]["chirality_status"] = inferred("unresolved", ev, why)
        m["stereochemistry"]["geometry_status"] = inferred("partially_specified", ev, why)
        m["stereochemistry"]["representation_status"] = inferred("not_encoded", ev, "Drawn geometry retained textually; full structure not encoded.")
    molecules[full] = m
    return full


def generic_mol(mid, label, page, fixed, variables):
    full = mol(mid, label, page, name=label)
    m = molecules[full]
    ev = [EV[page]]
    m["representation_type"] = "generic_template"
    m["smiles"] = F()
    m["generic_template"] = {
        "fixed_structure": F(fixed, {"description": fixed}, ev),
        "variables": [{
            "label": v,
            "attachment_points": F(f"source drawing attachment for {v}", {"label": v}, ev),
            "definition": F(desc, desc, ev),
            "allowed_values": F(allowed, allowed, ev),
            "constraints": F(),
        } for v, desc, allowed in variables],
        "instance_mapping": [],
        "reason": F("Generic R/Ar drawing in Scheme 1c", "The source explicitly uses variable substituents; no single concrete SMILES is valid.", ev),
    }
    m["stereochemistry"]["representation_status"] = inferred("not_encoded", ev, "Generic variable structure cannot encode one unique stereochemical graph.")
    return full


def q(number, unit, raw=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if raw is not None:
        value["raw_text"] = raw
    if qualifier:
        value["qualifier"] = qualifier
    return value


def classification(ev, changes=None, coverage="partial", why="Net bond changes inferred from source structures; no atom-mapped reaction was asserted."):
    c = st._skeleton("classification", ev)
    c["mapped_reaction_smiles"] = F()
    changes = changes or []
    c["bond_changes"] = [derived(x, ev, why) for x in changes]
    labels = {}
    for x in changes:
        if x["change_type"] in ("formed", "hydrogen_gain", "order_increased"):
            lab = st.bond_label(x["atoms"][0]["element"], x["atoms"][1]["element"], x["after_order"])
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(x["change_id"])
    c["labels"] = derived(list(labels.values()), ev, why) if coverage not in ("unresolved", "not_applicable") else F(None, [] if coverage == "not_applicable" else None, ev, method="inference", derivation=why) if coverage == "not_applicable" else F()
    c["coverage"] = inferred(coverage, ev, why)
    c["limitations"] = [
        "No complete atom mapping was asserted from the compact main-text drawing.",
        "Unprinted coproducts and proton sources remain unresolved.",
    ]
    c["model_votes"] = []
    return c


def formed(cid, a, b, order=1, count=None):
    x = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
         "before_order": 0, "after_order": order,
         "change_type": "hydrogen_gain" if b == "H" or a == "H" else "formed"}
    if count is not None:
        x["hydrogen_count"] = count
    return x


CO = [formed("b1", "C", "O")]


def rxn(rid, page, kind="specific_experiment", categories=None, name=None, transformation=None,
        classification_obj=None, extra_evidence=()):
    ev = [EV[page], *extra_evidence]
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = kind
    r["experiment_categories"] = F(categories or [], categories or [], ev)
    r["location"] = F(f"PDF page {page}", {"page": page, "region": rid}, ev)
    r["reaction_name"] = F(name, name, ev) if name else F()
    r["transformation_type"] = inferred(transformation, ev, "Transformation type summarized from the source arrow and text.") if transformation else F()
    r["classification"] = classification_obj or classification(ev, coverage="unresolved", why="Source structure mapping is unresolved; no formation label asserted.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["Machine visual extraction requires expert source review."]}
    r["completeness"] = "partial"
    return r


def participant(r, pid, mid, role, ev, *, stage="reaction", mass=None, amount=None,
                volume=None, equiv=None, molpct=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, ev)
    for key, value in (("mass", mass), ("amount", amount), ("volume", volume),
                       ("equivalents", equiv), ("mol_percent", molpct)):
        if value is not None:
            raw = value.get("raw_text") if isinstance(value, dict) else str(value)
            p["charges"][key] = F(raw or str(value), value, ev)
    if (equiv is not None or molpct is not None) and reference:
        p["equivalence_reference"] = inferred(reference, ev, "Loading/equivalent denominator interpreted as the explicitly stated limiting phenol/reactant.")
    r["participants"].append(p)
    return p


def conditions(r, ev, solvent=None, volume=None, temperature=None, time=None, reactor=None):
    if solvent:
        r["solvent_conditions"]["composition"] = F(solvent, [solvent], ev)
    if volume:
        r["solvent_conditions"]["total_volume"] = F(volume.get("raw_text", str(volume)), volume, ev)
    if reactor:
        r["special_conditions"]["reactor"] = F(reactor, reactor, ev)
    s = st.new_step("reaction", ev)
    s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    s["operation"] = F("reaction arrow", "reaction", ev)
    if temperature is not None:
        s["temperature"] = F(str(temperature), temperature, ev)
    if time is not None:
        s["time"] = F(time.get("raw_text", str(time)) if isinstance(time, dict) else str(time), time, ev)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, ev, value, *, yield_type=None, basis=None, scope=None, role="desired product"):
    o = st.new_outcome(oid, kind, targets, ev)
    o["product_role"] = F(role, role, ev)
    raw = value.get("raw_text", str(value)) if isinstance(value, dict) else str(value)
    o["value"] = F(raw, value, ev)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, ev)
    if basis:
        o["measurement_basis"] = F(basis, basis, ev)
    if scope:
        o["scope"] = F(scope, scope, ev)
    r["outcomes"].append(o)
    return o


# Core structures and common agents.
m1a = mol("1a", "1a", 2, "4-methoxyphenol", "COc1ccc(O)cc1")
m2a = mol("2a", "2a", 1, "2,3-butadienyl methyl carbonate", "C=C=CCOC(=O)OC")
m3aa = mol("3aa", "3aa", 2, "4-methoxyphenyl 2,3-butadienyl ether", "C=C=CCOc1ccc(OC)cc1")
m1b = mol("1b", "1b", 2, "phenol", "Oc1ccccc1")
m3ba = mol("3ba", "3ba", 2, "phenyl 2,3-butadienyl ether", "C=C=CCOc1ccccc1")
m1g = mol("1g", "1g", 2, "2-bromophenol", "Oc1ccccc1Br")
m3ga = mol("3ga", "3ga", 2, "2-bromophenyl 2,3-butadienyl ether", "C=C=CCOc1ccccc1Br")

agents = {}
def agent(key, label, page=2):
    agents[key] = mol("agent_" + key, label, page, label)
    return agents[key]

for key, label in {
    "pd_allyl_chloride_dimer": "[PdCl(π-allyl)]2", "pd2dba3": "Pd2dba3",
    "pdpph3_4": "Pd(PPh3)4", "binap": "BINAP", "xantphos": "Xantphos",
    "dpephos": "DPEPhos", "biphep": "BIPHEP", "sphos": "SPhos",
    "toluene": "toluene", "mtbe": "MTBE", "et2o": "Et2O", "dioxane": "dioxane",
    "thf": "THF", "ch3no2": "CH3NO2 internal standard",
}.items():
    agent(key, label)

# Scheme 1c generic current-work reaction.
gphenol = generic_mol("generic_phenol", "Ar-OH (1)", 1, "phenol/naphthol scaffold Ar-OH", [("Ar", "aryl group", "phenyl/naphthyl and substituted examples")])
gcarbonate = generic_mol("generic_allenyl_carbonate", "2 (R = H, alkyl or aryl)", 1, "R-substituted 2,3-allenylic methyl carbonate", [("R", "allene terminal substituent", ["H", "alkyl", "aryl"])])
gproduct = generic_mol("generic_allenyl_aryl_ether", "3 (R = H, alkyl or aryl)", 1, "allenic aromatic ether product", [("Ar", "aryl group", "phenyl/naphthyl and substituted examples"), ("R", "allene terminal substituent", ["H", "alkyl", "aryl"])])
r = rxn("scheme1c_general", 1, kind="generic_scheme", name="palladium-catalyzed C–O bond formation", transformation="phenol allenylic etherification", classification_obj=classification([EV[1]], CO))
participant(r, "phenol", gphenol, "reactant", [EV[1]])
participant(r, "carbonate", gcarbonate, "reactant", [EV[1]])
participant(r, "pd2dba3", agents["pd2dba3"], "catalyst", [EV[1]], molpct=q(2.5, "mol%", "2.5 mol%"), reference=gphenol)
participant(r, "xantphos", agents["xantphos"], "ligand", [EV[1]], molpct=q(10, "mol%", "10 mol%"), reference=gphenol)
participant(r, "et2o", agents["et2o"], "solvent", [EV[1]])
participant(r, "product", gproduct, "product", [EV[1]], product_role="generic desired product")
conditions(r, [EV[1]], "Et2O", temperature="rt")
r["review"]["unresolved"].append("Generic R/Ar scheme is not a concrete experiment and has no reported time or yield on page 1.")
bundle["reactions"].append(r)

# Preliminary current-work failed survey summarized across Scheme 1c/page-2 prose.
iodoanisole = mol("prelim_4_iodoanisole", "4-iodoanisole", 1, "4-iodoanisole")
allenols = generic_mol("prelim_allenols", "allenols", 1, "allenol(s)", [("R", "unspecified allenol substituent", "not supplied in the main text")])
prelim_prod = mol("prelim_desired_ether", "desired allenyl aryl ether product", 1, "desired allenyl aryl ether product")
r = rxn("preliminary_aryl_halide_allenol_survey", 1, categories=["control"], name="preliminary Pd/Cu-catalyzed coupling survey", transformation="attempted aryl ether formation", extra_evidence=(EV[2],))
participant(r, "aryl_halide", iodoanisole, "reactant", [EV[1], EV[2]])
participant(r, "allenols", allenols, "reactant", [EV[1], EV[2]])
participant(r, "desired_product", prelim_prod, "product", [EV[1], EV[2]], product_role="desired product")
conditions(r, [EV[1], EV[2]])
outcome(r, "trace", "trace", [prelim_prod], [EV[1], EV[2]], "trace", scope="summarized preliminary Pd- and Cu-catalyzed surveys")
r["review"]["unresolved"].append("The main text summarizes multiple Pd/Cu surveys but gives no concrete allenol identity or conditions; absent ESI prevents row-level reconstruction.")
bundle["reactions"].append(r)

# Table 1: eleven optimization rows plus the separate 1 mmol isolated experiment in entry 9 footnote b.
table1 = [
    (1, "pd_allyl_chloride_dimer", "binap", "toluene", 39, 48),
    (2, "pd2dba3", "binap", "toluene", 68, 26),
    (3, "pdpph3_4", "binap", "toluene", 43, 49),
    (4, "pd2dba3", "xantphos", "toluene", 80, 10),
    (5, "pd2dba3", "dpephos", "toluene", 28, 70),
    (6, "pd2dba3", "biphep", "toluene", 10, 89),
    (7, "pd2dba3", "sphos", "toluene", 68, 13),
    (8, "pd2dba3", "xantphos", "mtbe", 96, None),
    (9, "pd2dba3", "xantphos", "et2o", 99, None),
    (10, "pd2dba3", "xantphos", "dioxane", 67, 28),
    (11, "pd2dba3", "xantphos", "thf", 49, 47),
]
for entry, cat, lig, solv, yld, rec in table1:
    ev = [EV[2]]
    r = rxn(f"table1_entry_{entry}", 2, categories=["optimization"], name="optimization of Pd-catalyzed phenol–allenylic carbonate coupling", transformation="C–O bond formation", classification_obj=classification(ev, CO))
    participant(r, "1a", m1a, "reactant", ev, amount=q(0.2, "mmol", "0.2 mmol"))
    participant(r, "2a", m2a, "reactant", ev, equiv=q(1.5, "equiv", "1.5 equiv."), reference=m1a)
    participant(r, "pd", agents[cat], "catalyst", ev, molpct=q(5, "mol%", "[Pd] (5 mol%)"), reference=m1a)
    participant(r, "ligand", agents[lig], "ligand", ev, molpct=q(10, "mol%", "10 mol%"), reference=m1a)
    participant(r, "solvent", agents[solv], "solvent", ev, volume=q(1, "mL", "1 mL"))
    participant(r, "3aa", m3aa, "product", ev, product_role="desired product")
    r["reaction_scale"] = F("1a, 0.2 mmol", {"reference_molecule_id": m1a, "amount": q(0.2, "mmol")}, ev)
    conditions(r, ev, solv if solv != "et2o" else "Et2O", q(1, "mL", "1 mL"), "rt", q(12, "h", "12 h"))
    outcome(r, "yield", "yield", [m3aa], ev, q(yld, "%", f"{yld}%"), yield_type="NMR yield", basis="1H NMR analysis using CH3NO2 as internal standard", scope=f"Table 1 entry {entry}")
    if rec is not None:
        outcome(r, "recovery", "recovery", [m1a], ev, q(rec, "%", f"{rec}%"), yield_type="recovery", basis="1H NMR analysis using CH3NO2 as internal standard", scope=f"Table 1 entry {entry}", role="recovered starting material")
    else:
        r["notes"].append(F("Recovery column contains an em dash (—).", "source table dash; no numeric recovery assigned", ev))
    r["notes"].append(F("Table header reports [Pd] (5 mol%); for dimeric Pd sources this is retained as Pd loading rather than silently converted to complex mol%.", "Pd-atom loading expression retained", ev))
    bundle["reactions"].append(r)

r = rxn("table1_entry_9_1mmol_isolated", 2, categories=["optimization", "scale_up"], name="1 mmol isolated-yield repeat of Table 1 entry 9", transformation="C–O bond formation", classification_obj=classification([EV[2]], CO))
participant(r, "1a", m1a, "reactant", [EV[2]], amount=q(1, "mmol", "1 mmol"))
participant(r, "2a", m2a, "reactant", [EV[2]], equiv=q(1.5, "equiv", "1.5 equiv."), reference=m1a)
participant(r, "pd", agents["pd2dba3"], "catalyst", [EV[2]], molpct=q(5, "mol%", "[Pd] (5 mol%)"), reference=m1a)
participant(r, "ligand", agents["xantphos"], "ligand", [EV[2]], molpct=q(10, "mol%", "10 mol%"), reference=m1a)
participant(r, "solvent", agents["et2o"], "solvent", [EV[2]])
participant(r, "3aa", m3aa, "product", [EV[2]], product_role="desired product")
r["reaction_scale"] = F("1 mmol of 1a", {"reference_molecule_id": m1a, "amount": q(1, "mmol")}, [EV[2]])
conditions(r, [EV[2]], "Et2O", temperature="rt", time=q(12, "h", "12 h"))
outcome(r, "yield", "yield", [m3aa], [EV[2]], q(95, "%", "95%"), yield_type="isolated yield", basis="isolated", scope="Table 1 entry 9 footnote b, 1 mmol repeat")
r["review"]["unresolved"].append("Footnote b changes the scale to 1 mmol but does not explicitly restate solvent volume; the 1 mL header value was not assigned to this repeat.")
bundle["reactions"].append(r)

# Table 2 scope: 21 separate arrows.
scope2 = [
    ("1b", "3ba", 84, None), ("1c", "3ca", 94, None), ("1d", "3da", 71, None),
    ("1e", "3ea", 84, None), ("1f", "3fa", 90, None), ("1g", "3ga", 86, None),
    ("1h", "3ha", 92, None), ("1i", "3ia", 90, None), ("1j", "3ja", 94, None),
    ("1k", "3ka", 87, None), ("1l", "3la", 88, "b"), ("1m", "3ma", 82, "c"),
    ("1n", "3na", 70, "d"), ("1o", "3oa", 84, None), ("1p", "3pa", 95, None),
    ("1q", "3qa", 91, "e"), ("1r", "3ra", 94, "e"), ("1s", "3sa", 95, "e"),
    ("1t", "3ta", 90, "e"), ("1u", "3ua", 93, None), ("1v", "3va", 94, "f"),
]
names1 = {
    "1c": "p-cresol", "1d": "2,6-dimethylphenol", "1e": "2,6-diphenylphenol",
    "1f": "ortho-allyl-substituted phenol", "1h": "3-bromophenol", "1i": "4-bromophenol",
    "1j": "4-fluorophenol", "1k": "4-chlorophenol", "1l": "4-iodophenol",
    "1m": "methyl 4-hydroxybenzoate", "1n": "4-hydroxybenzonitrile",
    "1o": "naphthol 1o", "1p": "naphthol 1p", "1q": "hydroquinone", "1r": "resorcinol",
    "1s": "catechol", "1t": "1,1'-binaphthol", "1u": "hydroxymethylphenol 1u", "1v": "estrone",
}
scope2_mids = {"1b": m1b, "1g": m1g}
for s, _, _, _ in scope2:
    if s not in scope2_mids:
        scope2_mids[s] = mol(s, s, 2, names1[s], stereo="central_specified" if s == "1v" else "unresolved")
product2_mids = {"3ba": m3ba, "3ga": m3ga}
for _, p, _, _ in scope2:
    if p not in product2_mids:
        product2_mids[p] = mol(p, p, 2, f"Table 2 product {p}", stereo="central_specified" if p == "3va" else "unresolved")

for s, p, yld, note in scope2:
    ev = [EV[2]]
    sm, pm = scope2_mids[s], product2_mids[p]
    r = rxn(f"table2_{p}", 2, categories=["substrate_scope"], name="phenol/naphthol substrate-scope coupling", transformation="C–O bond formation", classification_obj=classification(ev, CO))
    participant(r, s, sm, "reactant", ev, amount=q(1, "mmol", "1.0 mmol"))
    two_amount = 2.5 if note == "e" else 1.5
    participant(r, "2a", m2a, "reactant", ev, amount=q(two_amount, "mmol", f"{two_amount} mmol"))
    catalyst = "pdpph3_4" if note in ("c", "d") else "pd2dba3"
    cat_loading = 5 if note in ("c", "d") else 2.5
    participant(r, "catalyst", agents[catalyst], "catalyst", ev, molpct=q(cat_loading, "mol%", f"{cat_loading} mol%"), reference=sm)
    if note not in ("c", "d"):
        participant(r, "xantphos", agents["xantphos"], "ligand", ev, molpct=q(10, "mol%", "10 mol%"), reference=sm)
    solvent_key = "mtbe" if note in ("b", "d", "f") else "dioxane" if note == "c" else "et2o"
    participant(r, "solvent", agents[solvent_key], "solvent", ev, volume=q(5, "mL", "5 mL"))
    participant(r, p, pm, "product", ev, product_role="desired product")
    temp = q(35, "°C", "35 °C") if note in ("b", "c") else q(40, "°C", "40 °C") if note == "d" else q(60, "°C", "60 °C") if note == "f" else "rt"
    conditions(r, ev, "MTBE" if solvent_key == "mtbe" else "dioxane" if solvent_key == "dioxane" else "Et2O", q(5, "mL", "5 mL"), temp, reactor="sealed tube" if note == "f" else None)
    r["reaction_scale"] = F("1 (1.0 mmol)", {"reference_molecule_id": sm, "amount": q(1, "mmol")}, ev)
    outcome(r, "yield", "yield", [pm], ev, q(yld, "%", f"{yld}%"), yield_type="isolated yield", basis="isolated", scope=f"Table 2 product {p}")
    r["review"]["unresolved"].append("Reaction time and workup/purification are not reported in the supplied main-text page; SI is not associated.")
    bundle["reactions"].append(r)

# Table 3 scope: 15 separate arrows. These substituted allenes can be axially chiral; no assignment is printed.
scope3 = [
    ("1b", "2b", "3bb", 88), ("1p", "2c", "3pc", 95), ("1p", "2d", "3pd", 98),
    ("1o", "2e", "3oe", 96), ("1b", "2f", "3bf", 99), ("1p", "2g", "3pg", 96),
    ("1p", "2h", "3ph", 98), ("1p", "2i", "3pi", 88), ("1b", "2j", "3bj", 71),
    ("1p", "2k", "3pk", 92), ("1p", "2l", "3pl", 75), ("1b", "2m", "3bm", 80),
    ("1p", "2n", "3pn", 75), ("1b", "2o", "3bo", 92), ("1p", "2p", "3pp", 90),
]
for s, c, p, yld in scope3:
    sm = scope2_mids[s]
    cm = mol(c, c, 3, f"2,3-allenylic carbonate {c}", stereo="axial_unspecified")
    pm = mol(p, p, 3, f"Table 3 allenic aromatic ether {p}", stereo="axial_unspecified")
    ev = [EV[3]]
    r = rxn(f"table3_{p}", 3, categories=["substrate_scope"], name="2,3-allenylic carbonate substrate-scope coupling", transformation="C–O bond formation", classification_obj=classification(ev, CO))
    participant(r, s, sm, "reactant", ev, amount=q(1, "mmol", "1.0 mmol"))
    participant(r, c, cm, "reactant", ev, amount=q(1.5, "mmol", "1.5 mmol"))
    participant(r, "pd2dba3", agents["pd2dba3"], "catalyst", ev, molpct=q(2.5, "mol%", "2.5 mol%"), reference=sm)
    participant(r, "xantphos", agents["xantphos"], "ligand", ev, molpct=q(10, "mol%", "10 mol%"), reference=sm)
    participant(r, "et2o", agents["et2o"], "solvent", ev, volume=q(5, "mL", "5 mL"))
    participant(r, p, pm, "product", ev, product_role="desired product")
    r["reaction_scale"] = F("1 (1.0 mmol)", {"reference_molecule_id": sm, "amount": q(1, "mmol")}, ev)
    conditions(r, ev, "Et2O", q(5, "mL", "5 mL"), "rt")
    outcome(r, "yield", "yield", [pm], ev, q(yld, "%", f"{yld}%"), yield_type="isolated yield", basis="isolated", scope=f"Table 3 product {p}")
    r["review"]["unresolved"].append("No axial stereochemical assignment, reaction time, workup, or purification is printed on the supplied page.")
    bundle["reactions"].append(r)

# Scheme 2a scale-up and four derivatizations.
r = rxn("scheme2a_gram_scale_3ga", 3, categories=["scale_up"], name="gram-scale allenylic ether synthesis", transformation="C–O bond formation", classification_obj=classification([EV[3]], CO))
participant(r, "1g", m1g, "reactant", [EV[3]], amount=q(5, "mmol", "5 mmol"))
participant(r, "2a", m2a, "reactant", [EV[3]], equiv=q(1.25, "equiv", "1.25 equiv."), reference=m1g)
participant(r, "pd2dba3", agents["pd2dba3"], "catalyst", [EV[3]], molpct=q(1.25, "mol%", "1.25 mol%"), reference=m1g)
participant(r, "xantphos", agents["xantphos"], "ligand", [EV[3]], molpct=q(5, "mol%", "5 mol%"), reference=m1g)
participant(r, "et2o", agents["et2o"], "solvent", [EV[3]])
participant(r, "3ga", m3ga, "product", [EV[3]], mass=q(1.0465, "g", "1.0465 g"), product_role="desired product")
r["reaction_scale"] = F("1g, 5 mmol", {"reference_molecule_id": m1g, "amount": q(5, "mmol")}, [EV[3]])
conditions(r, [EV[3]], "Et2O", temperature="rt")
outcome(r, "amount", "product_amount", [m3ga], [EV[3]], q(1.0465, "g", "1.0465 g"), scope="Scheme 2a")
outcome(r, "yield", "yield", [m3ga], [EV[3]], q(93, "%", "93%"), yield_type="yield", scope="Scheme 2a")
r["review"]["unresolved"].append("Reaction time, solvent volume, and yield measurement basis are not printed.")
bundle["reactions"].append(r)

# Derivatization helper.
def deriv_record(rid, product_label, product_name, changes, reagents, solvent, yld, stereo=None, selectivity=None):
    pm = mol(product_label, product_label, 3, product_name, stereo=stereo or "unresolved")
    ev = [EV[3]]
    r = rxn(rid, 3, categories=["derivatization"], name="synthetic transformation of the allene unit", transformation=product_name, classification_obj=classification(ev, changes))
    participant(r, "3ga", m3ga, "reactant", ev)
    for i, (key, label, role, dose) in enumerate(reagents):
        mid = agents.get(key) or agent(key, label, 3)
        kw = {}
        if dose and dose[1] == "mol%": kw["molpct"] = q(dose[0], "mol%", f"{dose[0]} mol%")
        elif dose and dose[1] == "equiv": kw["equiv"] = q(dose[0], "equiv", f"{dose[0]} equiv.")
        participant(r, f"reagent_{i+1}", mid, role, ev, reference=m3ga, **kw)
    skey = solvent.lower().replace("/", "_").replace(" ", "_")
    smid = agents.get(skey) or agent("solvent_" + skey, solvent, 3)
    participant(r, "solvent", smid, "solvent", ev)
    participant(r, product_label, pm, "product", ev, product_role="desired product")
    conditions(r, ev, solvent, temperature="rt")
    outcome(r, "yield", "yield", [pm], ev, q(yld, "%", f"{yld}%"), yield_type="yield", scope=f"Scheme 2a product {product_label}")
    if selectivity:
        outcome(r, "selectivity", "selectivity", [pm], ev, selectivity, scope=f"Scheme 2a; isomer order {selectivity.split('=')[0].strip()}", role="isomer mixture")
    r["review"]["unresolved"].append("Reaction time and measurement bases for yield/selectivity are not printed; no method was inferred.")
    bundle["reactions"].append(r)
    return pm

m4 = deriv_record("scheme2a_3ga_to_4", "4", "selenohydroxylation/oxidation product 4",
                  [formed("b1", "C", "Se"), formed("b2", "C", "O", 2)],
                  [("fluoropyridinium_otf", "1-fluoropyridinium triflate", "reagent", None), ("diphenyldiselane", "(PhSe)2", "reagent", None)],
                  "MeCN/THF (v/v, 1:1)", 42, stereo="geometry_mixture", selectivity="Z/E = 94/6")
m5 = deriv_record("scheme2a_3ga_to_5", "5", "bromohydroxylation product 5",
                  [formed("b1", "C", "Br"), formed("b2", "C", "O"), formed("b3", "O", "H")],
                  [("nbs", "NBS", "reagent", (1.2, "equiv"))], "1,4-dioxane/H2O (v/v, 1:1)", 72)
bundle["reactions"][-1]["notes"].append(F("sole product", "The source states product 5 was the sole product.", [EV[3]]))
m6 = deriv_record("scheme2a_3ga_to_6", "6", "gold-catalyzed hydration product 6",
                  [formed("b1", "C", "O"), formed("b2", "O", "H"), formed("b3", "C", "H")],
                  [("ipraucL", "IPrAuCl", "catalyst", (5, "mol%")), ("agotf", "AgOTf", "reagent", (5, "mol%"))],
                  "1,4-dioxane/H2O (v/v, 2:1)", 70, stereo="geometry_mixture", selectivity="E/Z = 98/2")
m7 = deriv_record("scheme2a_3ga_to_7", "7", "diiodide product (Z)-7",
                  [formed("b1", "C", "I"), formed("b2", "C", "I")],
                  [("iodine", "I2", "reagent", (1.2, "equiv"))], "Et2O", 69, stereo="geometry_Z")

# Scheme 2b eucalyptene A.
cinnamate = mol("methyl_4_hydroxycinnamate", "methyl 4-hydroxycinnamate", 3, "methyl 4-hydroxycinnamate", stereo="geometry_drawn")
eucalyptene = mol("eucalyptene_A", "eucalyptene A", 3, "eucalyptene A", stereo="geometry_drawn")
r = rxn("scheme2b_eucalyptene_A", 3, categories=["other"], name="synthesis of eucalyptene A", transformation="C–O bond formation", classification_obj=classification([EV[3]], CO))
participant(r, "cinnamate", cinnamate, "reactant", [EV[3]], amount=q(1, "mmol", "1 mmol"))
participant(r, "2a", m2a, "reactant", [EV[3]], equiv=q(1.5, "equiv", "1.5 equiv."), reference=cinnamate)
participant(r, "pdpph3", agents["pdpph3_4"], "catalyst", [EV[3]], molpct=q(5, "mol%", "5 mol%"), reference=cinnamate)
participant(r, "et2o", agents["et2o"], "solvent", [EV[3]])
participant(r, "eucalyptene", eucalyptene, "product", [EV[3]], product_role="natural product")
conditions(r, [EV[3]], "Et2O", temperature="rt")
outcome(r, "yield", "yield", [eucalyptene], [EV[3]], q(90, "%", "90%"), yield_type="yield", scope="Scheme 2b")
r["review"]["unresolved"].append("Reaction time, solvent volume, and yield measurement basis are unreported in the supplied pages.")
bundle["reactions"].append(r)

# Scheme 2c, one record per reaction arrow.
glucose = mol("D_glucose", "D-glucose", 3, "D-glucose", stereo="central_specified")
m8 = mol("8", "8", 3, "penta-O-trimethylsilyl-D-glucopyranose 8", stereo="central_specified")
m9 = mol("9", "9", 3, "glucosyl iodide 9", stereo="central_specified")
phenol_ch2oh = mol("p_hydroxymethylphenol", "p-hydroxymethylphenol", 3, "p-hydroxymethylphenol")
m10 = mol("10", "10", 3, "p-hydroxymethylphenyl 2,3-butadienyl ether 10")
terricollene = mol("terricollene_A", "terricollene A", 3, "terricollene A", stereo="central_specified")
for key, label in {"tmscl": "TMSCl", "pyridine": "pyridine", "tmsi": "TMSI", "dcm": "DCM", "dtbp": "2,6-di-tert-butylpyridine", "meoh": "MeOH", "resin_oh": "resin(OH−)"}.items():
    if key not in agents: agent(key, label, 3)

r = rxn("scheme2c_glucose_to_8", 3, categories=["other"], name="silyl protection of D-glucose", transformation="O-silylation", classification_obj=classification([EV[3]], [formed(f"b{i}", "O", "Si") for i in range(1, 6)]))
participant(r, "glucose", glucose, "reactant", [EV[3]], mass=q(4.0, "g", "4.0 g"))
participant(r, "tmscl", agents["tmscl"], "reagent", [EV[3]])
participant(r, "pyridine", agents["pyridine"], "solvent", [EV[3]])
participant(r, "8", m8, "product", [EV[3]], product_role="intermediate")
conditions(r, [EV[3]], "pyridine", temperature="rt", time="overnight")
r["review"]["unresolved"].append("TMSCl/pyridine quantities, product 8 yield, workup, and purification are unreported.")
bundle["reactions"].append(r)

r = rxn("scheme2c_8_to_9", 3, categories=["other"], name="formation of glucosyl iodide 9", transformation="anomeric iodination", classification_obj=classification([EV[3]], [formed("b1", "C", "I")]))
participant(r, "8", m8, "reactant", [EV[3]])
participant(r, "tmsi", agents["tmsi"], "reagent", [EV[3]])
participant(r, "dcm", agents["dcm"], "solvent", [EV[3]])
participant(r, "9", m9, "product", [EV[3]], product_role="intermediate")
conditions(r, [EV[3]], "DCM", temperature="rt", time=q(30, "min", "30 min"))
r["review"]["unresolved"].append("TMSI/DCM quantities, product 9 yield, workup, and purification are unreported.")
bundle["reactions"].append(r)

r = rxn("scheme2c_make_10", 3, categories=["scale_up"], name="70 mmol synthesis of allene 10", transformation="C–O bond formation", classification_obj=classification([EV[3]], CO))
participant(r, "phenol", phenol_ch2oh, "reactant", [EV[3]], amount=q(70, "mmol", "70 mmol"))
participant(r, "2a", m2a, "reactant", [EV[3]], equiv=q(1.14, "equiv", "1.14 equiv."), reference=phenol_ch2oh)
participant(r, "pd2dba3", agents["pd2dba3"], "catalyst", [EV[3]], molpct=q(1, "mol%", "1 mol%"), reference=phenol_ch2oh)
participant(r, "xantphos", agents["xantphos"], "ligand", [EV[3]], molpct=q(4, "mol%", "4 mol%"), reference=phenol_ch2oh)
participant(r, "et2o", agents["et2o"], "solvent", [EV[3]])
participant(r, "10", m10, "product", [EV[3]], mass=q(7.15, "g", "7.15 g"), product_role="synthetic intermediate")
conditions(r, [EV[3]], "Et2O", temperature="rt")
outcome(r, "amount", "product_amount", [m10], [EV[3], EV[4]], q(7.15, "g", "7.15 g"), scope="Scheme 2c")
r["review"]["unresolved"].append("Reaction time, solvent volume, isolated yield percentage, workup, and purification are unreported.")
bundle["reactions"].append(r)

r = rxn("scheme2c_9_plus_10_to_terricollene", 3, categories=["other"], name="glycosylation and desilylation to terricollene A", transformation="glycosidic C–O bond formation followed by desilylation", classification_obj=classification([EV[3], EV[4]], [formed("b1", "C", "O"), formed("b2", "O", "H", count=4)]), extra_evidence=(EV[4],))
participant(r, "9", m9, "reactant", [EV[3], EV[4]])
participant(r, "10", m10, "reactant", [EV[3], EV[4]], mass=q(7.15, "g", "7.15 g"))
participant(r, "base", agents["dtbp"], "reagent", [EV[3], EV[4]])
participant(r, "dcm", agents["dcm"], "solvent", [EV[3], EV[4]])
participant(r, "meoh", agents["meoh"], "reagent", [EV[3], EV[4]], stage="workup")
participant(r, "resin", agents["resin_oh"], "workup_agent", [EV[3], EV[4]], stage="workup")
participant(r, "terricollene", terricollene, "product", [EV[3], EV[4]], mass=q(4.5, "g", "4.5 g"), product_role="natural product")
conditions(r, [EV[3], EV[4]], "dry DCM", temperature="rt", time=q(4, "h", "4 h"))
s2 = st.new_step("desilylation", [EV[3], EV[4]], stage="workup")
s2["participant_ids"] = ["meoh", "resin", "terricollene"]
s2["operation"] = F("then MeOH, resin(OH−)", "methanol/resin desilylation", [EV[3], EV[4]])
s2["temperature"] = F("rt", "rt", [EV[3], EV[4]])
r["steps"].append(s2)
outcome(r, "amount", "product_amount", [terricollene], [EV[3], EV[4]], q(4.5, "g", "4.5 g"), scope="terricollene A isolated from the sequence")
outcome(r, "overall_yield", "yield", [terricollene], [EV[3], EV[4]], q(58, "%", "58% overall yield based on Glc"), yield_type="overall yield", scope="overall from D-glucose through intermediates 8 and 9 to terricollene A; not a last-step yield")
outcome(r, "anomer", "selectivity", [terricollene], [EV[3], EV[4]], "sole anomer", scope="anomeric selectivity; absolute assignment/method not printed", role="sole anomer product")
r["notes"].append(F("58% overall yield based on Glc", "Overall sequence yield begins at D-glucose and must not be interpreted as the isolated yield of only the 9 + 10 step.", [EV[3], EV[4]]))
r["review"]["unresolved"] += ["Amounts of 9, 2,6-di-tert-butylpyridine, DCM, MeOH, and resin are unreported.", "Duration of the methanol/resin stage and measurement bases for overall yield/anomeric selectivity are unreported."]
bundle["reactions"].append(r)

# Proposed mechanism is explicitly not a performed experiment.
meoh = agents["meoh"]
co2 = mol("CO2", "CO2", 4, "carbon dioxide")
pd0 = mol("Pd0", "Pd(0)", 4, "Pd(0)")
intermediate = generic_mol("mechanism_int", "Int", 4, "α-methylene-π-allylic Pd intermediate", [("R", "allene substituent", "H, alkyl or aryl")])
r = rxn("scheme3_proposed_mechanism", 4, kind="proposed_mechanism", name="proposed Pd(0)/α-methylene-π-allyl mechanism", transformation="proposed catalytic cycle", classification_obj=classification([EV[4]], coverage="not_applicable", why="Proposed elementary arrows are not classified as observed experimental bond formations."))
for pid, mid, role in [("phenol", gphenol, "reactant"), ("carbonate", gcarbonate, "reactant"), ("pd0", pd0, "catalyst"), ("int", intermediate, "unknown"), ("product", gproduct, "product"), ("co2", co2, "product"), ("meoh", meoh, "product")]:
    participant(r, pid, mid, role, [EV[4]], product_role="proposed species" if role == "product" else None)
r["notes"].append(F("A plausible mechanism is proposed", "Author-proposed mechanism, not an experimental verification.", [EV[4]]))
r["review"]["unresolved"].append("Elementary atom mapping, proton transfers, and intermediate structure details are proposed rather than experimentally established.")
bundle["reactions"].append(r)

bundle["molecules"] = list(molecules.values())
bundle["source_context"] = [
    {"context_id": f"{DOC}:scheme1a_natural_products", "document_id": DOC, "kind": "background_examples", "evidence_ids": [EV[1]], "note": "Pseudoxylarenone B, eucalyptene A, xyloallenolide A, and terricollene A are contextual natural-product structures, not experiments reported in this work."},
    {"context_id": f"{DOC}:scheme1b_prior_ether_synthesis", "document_id": DOC, "kind": "literature_context", "evidence_ids": [EV[1]], "note": "Generic Pd/Cu aryl-halide ether synthesis is cited prior work, not a this-paper experiment."},
    {"context_id": f"{DOC}:derivatization_citations", "document_id": DOC, "kind": "literature_context", "evidence_ids": [EV[3], EVT[3]], "note": "References 14–18 provide precedent for allene transformations; Scheme 2a arrows are this-work applications, while cited literature itself is not re-counted as experiments."},
    {"context_id": f"{DOC}:data_availability", "document_id": DOC, "kind": "document_note", "evidence_ids": [EV[4]], "note": "The article states that ESI contains experimental and analytical details, but no SI file is explicitly associated in source_mapping.json."},
    {"context_id": f"{DOC}:page5_references_only", "document_id": DOC, "kind": "references", "evidence_ids": [EV[5]], "note": "Page 5 contains references only and no new experiment arrow or table row from this work."},
]

# Keep the initialized submission file as the required template path, now populated.
dump(RUN / "submitted.json", bundle)

# Raw visual inventory and stereo diagnostics are kept apart from source-derived records.
diagnostics = {}
for mid in (m1a, m2a, m3aa, m1b, m3ba, m1g, m3ga):
    diagnostics[mid] = st.stereo_inspect(st.value(molecules[mid]["smiles"]))
dump(RUN / "raw_model/stereo_diagnostics.json", diagnostics)

inventory = {
    "document_id": DOC,
    "source": "visual inspection of supplied original page images",
    "pages": {
        "1": {"regions": ["title/abstract", "Scheme 1a natural products", "Scheme 1b cited prior work", "Scheme 1c preliminary failed approach and this-work generic method"], "reaction_boundaries": ["scheme1c_general", "preliminary_aryl_halide_allenol_survey"]},
        "2": {"regions": ["Table 1 rows 1-11 and footnotes", "Table 2 products 3ba-3va and footnotes", "discussion prose"], "reaction_boundaries": [f"table1_entry_{i}" for i in range(1, 12)] + ["table1_entry_9_1mmol_isolated"] + [f"table2_{p}" for _, p, _, _ in scope2]},
        "3": {"regions": ["Table 3 products 3bb-3pp and footnote", "Scheme 2a four transformations plus scale-up", "Scheme 2b eucalyptene A", "Scheme 2c four arrows"], "reaction_boundaries": [f"table3_{p}" for _, _, p, _ in scope3] + ["scheme2a_gram_scale_3ga", "scheme2a_3ga_to_4", "scheme2a_3ga_to_5", "scheme2a_3ga_to_6", "scheme2a_3ga_to_7", "scheme2b_eucalyptene_A", "scheme2c_glucose_to_8", "scheme2c_8_to_9", "scheme2c_make_10", "scheme2c_9_plus_10_to_terricollene"]},
        "4": {"regions": ["Scheme 3 proposed mechanism", "Scheme 2c prose continuation", "conclusions/data availability"], "reaction_boundaries": ["scheme3_proposed_mechanism"], "duplicate_evidence": ["scheme2c_make_10", "scheme2c_9_plus_10_to_terricollene"]},
        "5": {"regions": ["references 7-20"], "reaction_boundaries": []},
    },
    "limitations": ["No explicitly associated SI was supplied.", "Compact scope drawings not converted to full SMILES remain unresolved and needs_review.", "Page enumeration is not a claim of chemical accuracy."],
}
dump(RUN / "raw_model/page1_5_visual_inventory.json", inventory)

# Update task states without implying semantic accuracy.
tasks = st.read(RUN / "tasks.json")
for task in tasks:
    task["state"] = "attempted"
    task["reason"] = "Original page image inspected; detailed region-to-record mapping is in coverage/pages_0001_0005.json."
dump(RUN / "tasks.json", tasks)

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "source_context": len(bundle["source_context"])}, ensure_ascii=False))
