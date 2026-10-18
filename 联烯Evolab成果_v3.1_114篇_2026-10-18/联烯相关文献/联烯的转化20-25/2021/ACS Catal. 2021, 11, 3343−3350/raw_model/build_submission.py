#!/usr/bin/env python3
"""Author the bounded pages 1–6 v3.1 extraction from visually transcribed source facts."""
from pathlib import Path
import copy
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_a52c76191d5b058b"
RUN = ROOT / "runs/current"


def ev(page):
    return [f"{DOC}:page{page}:image", f"{DOC}:page{page}:text"]


def F(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None):
    if normalized is None and raw is not None:
        normalized = raw
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation)


def Q(number, unit, raw=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if raw is not None:
        value["raw_text"] = raw
    if qualifier:
        value["qualifier"] = qualifier
    return value


def molecule(mid, label, page, name=None, stereo_kind="unresolved", note=None,
             smiles=None, canonical=None, key=None, center=None):
    evidence = ev(page)
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(name or label, name or label, evidence)
    m["names"]["english"] = F(name or label, name or label, evidence)
    m["names"]["abbreviations"] = F()
    m["names"]["aliases"] = F()
    m["cas_numbers"] = F()
    m["identity_note"] = F(note, note, evidence) if note else F()
    m["external_properties"] = []
    deriv = "Stereochemical status assigned conservatively from the source drawing/label; absence of a descriptor was not treated as achiral or racemic."
    for field in ("atom_centers", "double_bonds", "axial", "mixture", "assignment"):
        m["stereochemistry"][field] = F()
    if stereo_kind == "product_branched":
        m["stereochemistry"]["chirality_status"] = F(None, "mixture", evidence, "inference", deriv)
        m["stereochemistry"]["geometry_status"] = F(None, "not_applicable", evidence, "inference", "The displayed branched products contain a terminal vinyl group for which E/Z is not applicable.")
        m["stereochemistry"]["representation_status"] = F(None, "not_encoded", evidence, "inference", "Wedge/dash information and reported dr/er were not converted to a guessed absolute stereochemical SMILES.")
        m["stereochemistry"]["mixture"] = F("reported dr and er values", "diastereomeric/enantiomeric mixture; ratios retained as experiment outcomes", evidence)
        m["stereochemistry"]["assignment"] = F("major relative stereoisomer drawn", "relative stereochemical drawing present; absolute configuration not assigned on pages 1–6", evidence)
    elif stereo_kind == "product_linear":
        m["stereochemistry"]["chirality_status"] = F(None, "unresolved", evidence, "inference", deriv)
        m["stereochemistry"]["geometry_status"] = F(None, "mixture", evidence, "inference", "The product label encompasses the reported E/Z ratio; no single-isomer identity was invented.")
        m["stereochemistry"]["representation_status"] = F(None, "not_encoded", evidence, "inference", "Source drawing retained by label; full stereochemical graph was not encoded.")
        m["stereochemistry"]["mixture"] = F("reported E/Z ratio", "E/Z mixture unless source states otherwise", evidence)
    elif stereo_kind == "specified_drawing":
        m["stereochemistry"]["chirality_status"] = F(None, "partially_specified", evidence, "inference", deriv)
        m["stereochemistry"]["geometry_status"] = F(None, "not_applicable", evidence, "inference", "No configurable alkene geometry in the displayed product identity.")
        m["stereochemistry"]["representation_status"] = F(None, "not_encoded", evidence, "inference", "Source wedge/dash drawing retained without guessed R/S encoding.")
        m["stereochemistry"]["assignment"] = F("wedge/dash structure shown", "relative/absolute drawing preserved; R/S descriptor not printed on pages 1–6", evidence)
    else:
        m["stereochemistry"]["chirality_status"] = F(None, "unresolved", evidence, "inference", deriv)
        m["stereochemistry"]["geometry_status"] = F(None, "unresolved", evidence, "inference", deriv)
        m["stereochemistry"]["representation_status"] = F(None, "unresolved", evidence, "inference", "No complete stereochemical graph was encoded from this label/drawing.")
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = F(f"source drawing of {label}", smiles, evidence, "structure_derivation",
                          "Candidate connectivity was manually transcribed from the page drawing; RDKit parsing/stereo diagnostics do not establish source correctness.")
        m["stereochemistry"]["isomeric_smiles"] = F(f"candidate for {label}", canonical or smiles, evidence, "structure_derivation",
                                                     "Canonical isomeric SMILES from standard stereo diagnostics.")
        m["stereochemistry"]["identity_key"] = F(None, key, evidence, "structure_derivation",
                                                   "stereo_graph_key generated by standard stereo diagnostics.")
        if center:
            m["stereochemistry"]["atom_centers"] = F(label, center, evidence, "structure_derivation",
                                                       "Center assignment comes from the encoded candidate graph; source label is retained separately.")
            m["stereochemistry"]["chirality_status"] = F(None, "specified", evidence, "structure_derivation",
                                                           "The source explicitly distinguishes L-1a and D-1a; the candidate graphs encode opposite centers.")
            m["stereochemistry"]["representation_status"] = F(None, "encoded", evidence, "structure_derivation",
                                                                "Candidate SMILES encodes the distinguished center.")
            m["stereochemistry"]["geometry_status"] = F(None, "not_applicable", evidence, "structure_derivation",
                                                          "Only a carbonyl double bond is present; E/Z does not apply.")
        elif mid == "m_1a":
            m["stereochemistry"]["atom_centers"] = F("stereocenter drawn without an R/S label in the scope", [{"atom_index": 3, "assignment": "unspecified"}], evidence)
            m["stereochemistry"]["chirality_status"] = F(None, "unspecified", evidence, "structure_derivation",
                                                           "The connectivity has one potential center, but pages 1–6 do not assign R/S to the generic 1a identity.")
            m["stereochemistry"]["geometry_status"] = F(None, "not_applicable", evidence, "structure_derivation", "No E/Z-capable alkene is present.")
            m["stereochemistry"]["representation_status"] = F(None, "partial", evidence, "structure_derivation", "Connectivity is encoded while the potential center remains unspecified.")
        elif mid in {"m_2a", "m_terminal_allene_oh"}:
            m["stereochemistry"]["chirality_status"] = F(None, "achiral", evidence, "structure_derivation", "The terminal CH2 allene candidate cannot possess configurationally stable allene axial chirality.")
            m["stereochemistry"]["geometry_status"] = F(None, "not_applicable", evidence, "structure_derivation", "Cumulated terminal allene geometry is not represented as E/Z.")
            m["stereochemistry"]["representation_status"] = F(None, "partial", evidence, "structure_derivation", "Connectivity is encoded; toolkit allene/axial reporting is limited.")
            m["stereochemistry"]["axial"] = F("terminal allene drawing", "not axially chiral because one terminus is CH2", evidence, "structure_derivation", "Derived from terminal substitution pattern.")
        elif mid == "m_ethyl_acrylate":
            m["stereochemistry"]["chirality_status"] = F(None, "achiral", evidence, "structure_derivation", "Candidate graph contains no stereogenic center.")
            m["stereochemistry"]["geometry_status"] = F(None, "not_applicable", evidence, "structure_derivation", "Terminal alkene does not have E/Z geometry.")
            m["stereochemistry"]["representation_status"] = F(None, "encoded", evidence, "structure_derivation", "Candidate graph is fully encoded for applicable stereochemistry.")
    return m


MOLECULES = {}


def ensure(mid, label, page, **kwargs):
    if mid not in MOLECULES:
        MOLECULES[mid] = molecule(mid, label, page, **kwargs)
    else:
        occ = MOLECULES[mid]["occurrences"]
        if not any(o["evidence_ids"] == ev(page) for o in occ):
            occ.append({"document_id": DOC, "compound_label": F(label, label, ev(page)), "evidence_ids": ev(page)})
    return mid


def part(pid, mid, role, evidence, ref=None, amount=None, equiv=None, molpct=None, volume=None,
         stage="reaction", product_role=None, inherited=False):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, evidence)
    method = "inherited_general_condition" if inherited else "original_extraction"
    if amount:
        p["charges"]["amount"] = F(amount[0], amount[1], evidence, method)
    if equiv:
        p["charges"]["equivalents"] = F(equiv[0], equiv[1], evidence, method)
    if molpct:
        p["charges"]["mol_percent"] = F(molpct[0], molpct[1], evidence, method)
    if volume:
        p["charges"]["volume"] = F(volume[0], volume[1], evidence, method)
    if ref:
        p["equivalence_reference"] = F(ref, ref, evidence, method)
    return p


def class_partial(evidence, labels):
    changes = []
    expected = []
    for i, code in enumerate(labels, 1):
        if code == "C-C": atoms = [{"element": "C", "atom_map": None}, {"element": "C", "atom_map": None}]
        elif code == "C-H": atoms = [{"element": "C", "atom_map": None}, {"element": "H", "atom_map": None}]
        elif code == "C-O": atoms = [{"element": "C", "atom_map": None}, {"element": "O", "atom_map": None}]
        elif code == "O-H": atoms = [{"element": "O", "atom_map": None}, {"element": "H", "atom_map": None}]
        elif code == "C-N": atoms = [{"element": "C", "atom_map": None}, {"element": "N", "atom_map": None}]
        else: raise ValueError(code)
        cid = f"b{i}"
        ch = {"change_id": cid, "atoms": atoms, "before_order": 0, "after_order": 1,
              "change_type": "hydrogen_gain" if "H" in [a["element"] for a in atoms] else "formed"}
        changes.append(F(None, ch, evidence, "inference", "Net bond formation inferred from the displayed reactant/product drawings; atom mapping is not claimed."))
        lab = st.bond_label(atoms[0]["element"], atoms[1]["element"], 1)
        lab["supporting_change_ids"] = [cid]
        expected.append(lab)
    return {"standard_id": st.STANDARD, "mapped_reaction_smiles": F(), "bond_changes": changes,
            "labels": F(None, expected, evidence, "inference", "Taxonomy labels derived from the partial net bond-change description."),
            "coverage": F(None, "partial", evidence, "inference", "Full atom mapping, coproduct balance, and proton source are unresolved."),
            "limitations": ["Partial net classification from source drawings; no complete atom mapping.", "CO2 and other coproducts are not represented."],
            "model_votes": []}


def class_na(evidence, reason):
    return {"standard_id": st.STANDARD, "mapped_reaction_smiles": F(), "bond_changes": [],
            "labels": F(None, [], evidence, "inference", reason),
            "coverage": F(None, "not_applicable", evidence, "inference", reason),
            "limitations": [reason], "model_votes": []}


def class_unresolved(evidence, reason):
    return {"standard_id": st.STANDARD, "mapped_reaction_smiles": F(), "bond_changes": [],
            "labels": F(), "coverage": F(None, "unresolved", evidence, "inference", reason),
            "limitations": [reason], "model_votes": []}


def base_reaction(rid, page, rtype, category, location, name, transform, evidence=None):
    evidence = evidence or ev(page)
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = rtype
    r["experiment_categories"] = F(category, category, evidence)
    r["location"] = F(location, {"page": page, "region": location}, evidence)
    r["reaction_name"] = F(name, name, evidence) if name else F()
    r["transformation_type"] = F(None, transform, evidence, "inference", f"Transformation type inferred from source caption/drawing: {transform}.") if transform else F()
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["Supporting Information was mentioned but not supplied or associated; omitted experimental details remain unreported."]}
    r["completeness"] = "partial"
    return r


def reaction_step(r, evidence, pids, operation, time=None, temperature=None):
    s = st.new_step("reaction", evidence)
    s["participant_ids"] = list(pids)
    s["operation"] = F(operation, operation, evidence)
    if time:
        s["time"] = F(time[0], time[1], evidence, "inherited_general_condition" if time[2] else "original_extraction")
    if temperature:
        s["temperature"] = F(temperature[0], temperature[1], evidence, "inherited_general_condition" if temperature[2] else "original_extraction")
    r["steps"] = [s]


def outcome(r, oid, kind, targets, raw, normalized, evidence, basis_raw=None, scope=None, yield_type=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = F(raw, normalized, evidence)
    if basis_raw:
        o["measurement_basis"] = F(basis_raw, basis_raw, evidence)
    if scope:
        o["scope"] = F(scope, scope, evidence)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, evidence)
    o["product_role"] = F("main product", "main product", evidence) if targets else F()
    r["outcomes"].append(o)


# Core identities and shared agents.
ensure("m_1a", "1a", 3, name="N-phenyl-protected phenylalanine",
       smiles="O=C(O)C(Nc1ccccc1)Cc1ccccc1", canonical="O=C(O)C(Cc1ccccc1)Nc1ccccc1",
       key="c4ef5cccd9be2733a11eca2f428fd8b42a2543b557029272198119c9c9ea8b6e")
ensure("m_L1a", "L-1a", 4, name="N-phenyl amino acid L-1a",
       smiles="O=C(O)[C@H](Nc1ccccc1)Cc1ccccc1", canonical="O=C(O)[C@@H](Cc1ccccc1)Nc1ccccc1",
       key="37de646ef2b7922b2b26018028eef8cf13d0c2717d007bdb18b49831906287ca",
       center=[{"atom_index": 3, "configuration": "R in candidate graph"}])
ensure("m_D1a", "D-1a", 4, name="N-phenyl amino acid D-1a",
       smiles="O=C(O)[C@@H](Nc1ccccc1)Cc1ccccc1", canonical="O=C(O)[C@H](Cc1ccccc1)Nc1ccccc1",
       key="9f9fe4baa19ca56e2fc821925d121269aa383ca8a2ef2257a7e7d1b30cbb25d0",
       center=[{"atom_index": 3, "configuration": "S in candidate graph"}])
ensure("m_2a", "2a", 3, name="cyclohexyloxyallene", smiles="C=C=COC1CCCCC1",
       canonical="C=C=COC1CCCCC1", key="51f4098b436bbb176bb951dfabdd3965ba5789551c059584deb5b869253dde92")
ensure("m_terminal_allene_oh", "terminal allene bearing (CH2)3OH", 6,
       smiles="C=C=CCCCO", canonical="C=C=CCCCO", key="4ac596f91feec7ef7c9616be16f4502d0e8db5d3f6540ff07939a5e986e349e9")
ensure("m_ethyl_acrylate", "ethyl acrylate", 4, name="ethyl acrylate", smiles="C=CC(=O)OCC",
       canonical="C=CC(=O)OCC", key="141f9b79545c089753b37ff4a7296df306000a690df71dbe841e7d32ccd3da9e")

for mid, label, page, name in [
    ("m_pd2dba3", "[Pd2(dba)3]", 3, "tris(dibenzylideneacetone)dipalladium"),
    ("m_L9", "L9", 3, "Feringa ligand L9"), ("m_entL9", "ent-L9", 5, "enantiomeric Feringa ligand ent-L9"),
    ("m_sbinap", "(S)-BINAP", 3, "(S)-BINAP"), ("m_racbinap", "rac-BINAP", 6, "rac-BINAP"),
    ("m_acid", "p-CF3C6H4COOH", 3, "4-(trifluoromethyl)benzoic acid"),
    ("m_pc1", "PC1", 3, "[Ir(ppy)2(dtbbpy)]PF6"),
    ("m_pc2", "PC2", 3, "[Ir(dF(CF3)ppy)2(dtbbpy)]PF6"),
    ("m_pc3", "PC3", 3, "4-CzIPN"), ("m_dmf", "DMF", 3, "N,N-dimethylformamide"),
    ("m_mecn", "MeCN", 6, "acetonitrile"), ("m_tempo", "TEMPO", 4, "TEMPO"),
    ("m_cs2co3", "Cs2CO3", 5, "cesium carbonate"), ("m_cyclohexyl_allene", "cyclohexyl allene", 4, "cyclohexyl allene"),
    ("m_9bbn", "9-BBN", 6, "9-BBN"), ("m_h2o2", "H2O2", 6, "hydrogen peroxide"),
    ("m_tfa", "TFA", 6, "trifluoroacetic acid"), ("m_dcm", "DCM", 6, "dichloromethane"),
    ("m_pdc", "Pd/C", 6, "palladium on carbon"), ("m_h2", "H2", 6, "hydrogen"),
    ("m_meoh", "MeOH", 6, "methanol"), ("m_etoac", "EtOAc", 6, "ethyl acetate"),
    ("m_rh", "[Rh(acac)(CO)2]", 6, "dicarbonylacetylacetonato rhodium(I)"),
    ("m_biphephos", "BIPHEPHOS", 6, "BIPHEPHOS"), ("m_co", "CO", 6, "carbon monoxide"),
    ("m_thf", "THF", 6, "tetrahydrofuran")]:
    ensure(mid, label, page, name=name)

REACTIONS = []

# Page 2 conceptual and prior-work schemes.
for rid, rtype, loc, title in [
    ("scheme1a_retrosynthesis", "generic_scheme", "Scheme 1A", "Retrosynthetic analysis of 1,2-amino alcohols"),
    ("scheme1b_prior_alkoxyallene", "literature_example", "Scheme 1B", "Reported soft carbon nucleophile additions to alkoxyallenes"),
    ("scheme1c_this_work", "generic_scheme", "Scheme 1C", "Regioselective decarboxylative hydroaminoalkylation of allenes")]:
    r = base_reaction(rid, 2, rtype, [], loc, title, None)
    r["classification"] = class_unresolved(ev(2), "Generic variables and/or cited prior-work structures prevent complete reaction mapping.")
    r["review"]["unresolved"].append("Conceptual/generic scheme; not counted as a concrete experiment performed under a fully specified condition set.")
    REACTIONS.append(r)


def branched_participants(r, evidence, aa_mid, allene_mid, pc_mid="m_pc2", ligand_mid="m_L9",
                          pc_loading=1.0, include_pd=True, include_ligand=True, inherited=False):
    condition_evidence = list(dict.fromkeys(evidence + (ev(3) if inherited else [])))
    aa = part("aa", aa_mid, "reactant", condition_evidence, ref=aa_mid,
              amount=("0.1 mmol", Q(0.1, "mmol")), equiv=("1 equiv", Q(1, "equiv")), inherited=inherited)
    al = part("allene", allene_mid, "reactant", condition_evidence, ref=aa_mid,
              amount=("0.15 mmol", Q(0.15, "mmol")), equiv=("1.5 equiv", Q(1.5, "equiv")), inherited=inherited)
    r["participants"] += [aa, al]
    if include_pd:
        r["participants"].append(part("pd", "m_pd2dba3", "catalyst", condition_evidence, ref=aa_mid,
                                      molpct=("2.5 mol %", Q(2.5, "mol%")), inherited=inherited))
    if include_ligand:
        r["participants"].append(part("ligand", ligand_mid, "ligand", condition_evidence, ref=aa_mid,
                                      molpct=(("6 mol %" if ligand_mid == "m_sbinap" else "12 mol %"), Q(6 if ligand_mid == "m_sbinap" else 12, "mol%")), inherited=inherited))
    r["participants"].append(part("acid", "m_acid", "additive", condition_evidence, ref=aa_mid,
                                  molpct=("50 mol %", Q(50, "mol%")), inherited=inherited))
    if pc_mid:
        r["participants"].append(part("pc", pc_mid, "catalyst", condition_evidence, ref=aa_mid,
                                      molpct=(f"{pc_loading:g} mol %", Q(pc_loading, "mol%")), inherited=inherited))
    r["participants"].append(part("solvent", "m_dmf", "solvent", condition_evidence,
                                  volume=("1 mL", Q(1, "mL")), inherited=inherited))
    r["reaction_scale"] = F("amino acid, 0.1 mmol", {"reference_molecule_id": aa_mid, "amount": Q(0.1, "mmol")}, condition_evidence,
                            "inherited_general_condition" if inherited else "original_extraction")
    r["solvent_conditions"]["composition"] = F("DMF", ["DMF"], condition_evidence, "inherited_general_condition" if inherited else "original_extraction")
    r["solvent_conditions"]["total_volume"] = F("1 mL", Q(1, "mL"), condition_evidence, "inherited_general_condition" if inherited else "original_extraction")
    r["solvent_conditions"]["concentration"] = F("0.1 M", Q(0.1, "M"), condition_evidence, "inherited_general_condition" if inherited else "original_extraction")
    r["solvent_conditions"]["concentration_basis"] = F("Table header/source condition", "nominal reaction concentration", condition_evidence, "inherited_general_condition" if inherited else "original_extraction")


def add_branched_conditions(r, evidence, pids, light=True, inherited=False):
    condition_evidence = list(dict.fromkeys(evidence + (ev(3) if inherited else [])))
    reaction_step(r, condition_evidence, pids, "irradiation under the stated dual Pd/photoredox conditions" if light else "reaction without blue LED",
                  ("16 h", Q(16, "h"), inherited), ("room temperature", "room temperature", inherited))
    if light:
        r["special_conditions"]["light_source"] = F("blue LED", "blue LED", condition_evidence, "inherited_general_condition" if inherited else "original_extraction")


# Initial (S)-BINAP experiment described in prose before Table 1.
r = base_reaction("initial_sbinap_3aa", 3, "specific_experiment", ["optimization"], "page 3 prose above Table 1",
                  "decarboxylative hydroaminoalkylation", "branched decarboxylative hydroaminoalkylation")
ensure("m_3aa", "3aa", 3, stereo_kind="product_branched")
branched_participants(r, ev(3), "m_1a", "m_2a", "m_pc1", "m_sbinap", 1)
r["participants"].append(part("product", "m_3aa", "product", ev(3), product_role="main product"))
add_branched_conditions(r, ev(3), [p["participant_id"] for p in r["participants"]])
outcome(r, "yield", "yield", ["m_3aa"], "good yields", "good yields", ev(3), "basis not stated in prose", yield_type="reported qualitative yield")
for oid, raw, normalized, scope in [("bl", ">19:1 B/L", ">19:1", "branched/linear ratio"), ("dr", "1.3:1 dr", "1.3:1", "diastereomeric ratio"), ("er", "69:31 er", "69:31", "enantiomeric ratio")]:
    outcome(r, oid, "selectivity", ["m_3aa"], raw, normalized, ev(3), "dr by crude 1H NMR; er by chiral HPLC" if oid != "bl" else None, scope=scope)
r["classification"] = class_partial(ev(3), ["C-C", "C-H"])
REACTIONS.append(r)


opt_rows = [
    (1, "m_pc1", 1, 72, None, "86:14", "4:1", "97:3/89:11", None),
    (2, "m_pc2", 1, 73, 70, "90:10", "4.5:1", "96:4/88:12", None),
    (3, "m_pc3", 1, 41, None, "85:15", "3.2:1", "96:4/88:12", None),
    (4, "m_pc2", 2, 67, None, "89:11", "4.5:1", "96:4/90:10", None),
    (5, "m_pc2", 0.5, 56, None, "90:10", "4.2:1", "96:4/88:12", None),
    (6, None, None, None, None, None, None, None, "without photocatalyst"),
    (7, "m_pc2", 1, None, None, None, None, None, "without blue LED"),
    (8, "m_pc2", 1, None, None, None, None, None, "without Pd and ligand"),
]
for entry, pc, loading, nmr_yield, iso_yield, bl, dr, eratio, modification in opt_rows:
    evidence = ev(3) + (ev(5) if entry == 2 else [])
    r = base_reaction(f"table1_entry{entry}", 3, "specific_experiment", ["optimization", "control"] if entry >= 6 else ["optimization"],
                      f"Table 1 entry {entry}", "decarboxylative hydroaminoalkylation", "branched decarboxylative hydroaminoalkylation", evidence)
    ensure("m_3aa", "3aa", 5 if entry == 2 else 3, stereo_kind="product_branched")
    branched_participants(r, evidence, "m_1a", "m_2a", pc, "m_L9", loading or 0,
                          include_pd=entry != 8, include_ligand=entry != 8)
    r["participants"].append(part("product", "m_3aa", "product", evidence, product_role="main product"))
    add_branched_conditions(r, evidence, [p["participant_id"] for p in r["participants"]], light=entry != 7)
    if modification:
        r["procedure"]["modifications"] = F(modification, modification, ev(3))
    if entry >= 6:
        outcome(r, "nr", "no_reaction", [], "nr", "no reaction", ev(3), "detection method not stated")
        r["classification"] = class_na(ev(3), "Negative control reported no reaction; no intended bond formation is asserted.")
    else:
        outcome(r, "nmr_yield", "yield", ["m_3aa"], f"{nmr_yield}%", Q(nmr_yield, "%"), ev(3), "NMR yield using isoquinoline as internal standard", yield_type="NMR yield")
        if iso_yield:
            outcome(r, "isolated_yield", "yield", ["m_3aa"], f"({iso_yield})", Q(iso_yield, "%"), evidence, "isolated yield (parentheses in Table 1; repeated as 70% in Table 3)", yield_type="isolated yield")
        for oid, raw, scope, basis in [("bl", bl, "branched/linear ratio", None), ("dr", dr, "diastereomeric ratio", "1H NMR of crude mixture"), ("er", eratio, "enantiomeric ratio major/minor diastereomer", "chiral HPLC")]:
            outcome(r, oid, "selectivity", ["m_3aa"], f"{raw} {oid if oid != 'bl' else 'B/L'}", raw, evidence, basis, scope=scope)
        r["classification"] = class_partial(evidence, ["C-C", "C-H"])
    REACTIONS.append(r)


table2 = [
    ("b",70,"3.8:1","87:13","96:4/89:11",False,False),("c",63,"5:1","91:9","96:4/90:10",True,False),
    ("d",74,"4.2:1","87:13","96:4/87:13",True,False),("e",59,"3.6:1","93:7","96:4/92:8",False,False),
    ("f",41,"3.1:1","90:10","94:6/86:14",True,False),("g",51,"3.6:1","87:13","96:4/88:12",True,True),
    ("h",49,"4:1","87:13","96:4/89:11",True,True),("i",45,"4:1","88:12","95:5/74:26",True,False),
    ("j",55,"4.3:1","88:12","97:3/90:10",True,True),("k",54,"3.2:1","89:11","94:6/81:19",True,False),
    ("l",63,"4.5:1","87:13","95:5/81:19",False,False),("m",49,"3.5:1","88:12","96:4/87:13",False,False),
    ("n",58,"5.3:1","88:12","95:5/84:16",True,False),("o",53,"3.8:1","89:11","96:4/85:15",True,False),
    ("p",44,"4.1:1","84:16","93:7/85:15",False,False),
]
for letter, y, dr, bl, eratio, use_pc1, deriv_er in table2:
    aa_mid = ensure(f"m_1{letter}", f"1{letter}", 4, note="Amino-acid substrate identity retained from the Table 2 source drawing; no complete SMILES asserted.")
    prod_mid = ensure(f"m_3{letter}a", f"3{letter}a", 4, stereo_kind="product_branched")
    r = base_reaction(f"table2_3{letter}a", 4, "specific_experiment", ["substrate_scope"], f"Table 2 product 3{letter}a",
                      "decarboxylative hydroaminoalkylation", "branched decarboxylative hydroaminoalkylation")
    branched_participants(r, ev(4), aa_mid, "m_2a", "m_pc1" if use_pc1 else "m_pc2", inherited=True)
    r["participants"].append(part("product", prod_mid, "product", ev(4), product_role="main product"))
    add_branched_conditions(r, ev(4), [p["participant_id"] for p in r["participants"]], inherited=True)
    outcome(r, "yield", "yield", [prod_mid], f"{y}% yield", Q(y, "%"), ev(4), "isolated yield", yield_type="isolated yield")
    outcome(r, "dr", "selectivity", [prod_mid], f"{dr} dr", dr, ev(4), "1H NMR of crude mixture", scope="diastereomeric ratio")
    outcome(r, "bl", "selectivity", [prod_mid], f"{bl} B/L", bl, ev(4), None, scope="branched/linear ratio")
    erbasis = "chiral HPLC after derivatization" if deriv_er else "chiral HPLC"
    outcome(r, "er", "selectivity", [prod_mid], f"{eratio} er", eratio, ev(4), erbasis, scope="enantiomeric ratio major/minor diastereomer")
    r["classification"] = class_partial(ev(4), ["C-C", "C-H"])
    REACTIONS.append(r)


table3 = [
    ("b",78,"4.3:1","91:9","96:4/85:15",False),("c",71,"5.6:1","91:9","97:3/87:13",False),
    ("d",75,"3.4:1",">95:5","90:10/74:26",False),("e",73,"2.7:1",">95:5","89:11/74:26",False),
    ("f",81,"3.9:1","88:12","91:9/72:28",False),("g",81,"3.8:1","91:9","92:8/75:25",False),
    ("h",81,"4.5:1","87:13","90:10/69:31",False),("i",77,"4.7:1","88:12","92:8/76:24",False),
    ("j",53,"4.5:1","91:9","96:4/90:10",False),("k",57,"5.7:1","88:12","96:4/93:7",False),
    ("m",37,"2:1",">19:1","95:5/89:11",True),
]
for letter, y, dr, bl, eratio, use_pc1 in table3:
    allene_mid = ensure(f"m_2{letter}", f"2{letter}", 5, note="Alkoxyallene identity retained from the Table 3 source drawing; no complete SMILES asserted.")
    prod_mid = ensure(f"m_3a{letter}", f"3a{letter}", 5, stereo_kind="product_branched")
    r = base_reaction(f"table3_3a{letter}", 5, "specific_experiment", ["substrate_scope"], f"Table 3 product 3a{letter}",
                      "decarboxylative hydroaminoalkylation", "branched decarboxylative hydroaminoalkylation")
    branched_participants(r, ev(5), "m_1a", allene_mid, "m_pc1" if use_pc1 else "m_pc2", inherited=True)
    r["participants"].append(part("product", prod_mid, "product", ev(5), product_role="main product"))
    add_branched_conditions(r, ev(5), [p["participant_id"] for p in r["participants"]], inherited=True)
    outcome(r, "yield", "yield", [prod_mid], f"{y}% yield", Q(y, "%"), ev(5), "isolated yield", yield_type="isolated yield")
    outcome(r, "dr", "selectivity", [prod_mid], f"{dr} dr", dr, ev(5), "1H NMR of crude mixture", scope="diastereomeric ratio")
    outcome(r, "bl", "selectivity", [prod_mid], f"{bl} B/L", bl, ev(5), None, scope="branched/linear ratio")
    outcome(r, "er", "selectivity", [prod_mid], f"{eratio} er", eratio, ev(5), "chiral HPLC", scope="enantiomeric ratio major/minor diastereomer")
    r["classification"] = class_partial(ev(5), ["C-C", "C-H"])
    REACTIONS.append(r)

# Matched/mismatched chiral-alkoxyallene entries 3al and 3al'.
ensure("m_2l", "2l", 5, stereo_kind="specified_drawing", note="Chiral alkoxyallene source drawing retained; absolute descriptor not printed on page 5.")
for suffix, product_label, ligand, y, synanti, sync3, antic3 in [
    ("entL9", "3al", "m_entL9", 83, "7.3:1", ">19:1", "12:1"),
    ("L9", "3al'", "m_L9", 60, "7.1:1", "2.2:1", "2:1")]:
    pmid = ensure("m_3al" if suffix == "entL9" else "m_3al_prime", product_label, 5, stereo_kind="specified_drawing")
    r = base_reaction(f"table3_{product_label.replace(chr(39),'prime')}_{suffix}", 5, "specific_experiment", ["substrate_scope"], f"Table 3 {product_label} ({suffix})",
                      "decarboxylative hydroaminoalkylation", "branched decarboxylative hydroaminoalkylation")
    branched_participants(r, ev(5), "m_1a", "m_2l", "m_pc2", ligand, inherited=True)
    r["participants"].append(part("product", pmid, "product", ev(5), product_role="main product"))
    add_branched_conditions(r, ev(5), [p["participant_id"] for p in r["participants"]], inherited=True)
    outcome(r, "yield", "yield", [pmid], f"{y}% yield", Q(y, "%"), ev(5), "isolated yield", yield_type="isolated yield")
    for oid, raw, sc in [("synanti", synanti, "1,2-syn:anti dr"), ("sync3", sync3, "1,2-syn,C3 dr"), ("antic3", antic3, "1,2-anti,C3 dr"), ("bl", "91:9", "branched/linear ratio")]:
        outcome(r, oid, "selectivity", [pmid], f"{raw} {sc}", raw, ev(5), "1H NMR of crude mixture" if oid != "bl" else None, scope=sc)
    r["classification"] = class_partial(ev(5), ["C-C", "C-H"])
    REACTIONS.append(r)


def terminal_scope_participants(r, evidence, aa_mid, allene_mid):
    r["participants"] = [
        part("aa", aa_mid, "reactant", evidence, ref=aa_mid, amount=("0.2 mmol", Q(0.2, "mmol"))),
        part("allene", allene_mid, "reactant", evidence, ref=aa_mid, equiv=("1.2 equiv", Q(1.2, "equiv"))),
        part("pd", "m_pd2dba3", "catalyst", evidence, ref=aa_mid, molpct=("2.5 mol %", Q(2.5, "mol%"))),
        part("ligand", "m_racbinap", "ligand", evidence, ref=aa_mid, molpct=("6 mol %", Q(6, "mol%"))),
        part("pc", "m_pc1", "catalyst", evidence, ref=aa_mid, molpct=("1 mol %", Q(1, "mol%"))),
        part("solvent", "m_mecn", "solvent", evidence),
    ]
    r["reaction_scale"] = F("1 (0.2 mmol)", {"reference_molecule_id": aa_mid, "amount": Q(0.2, "mmol")}, evidence)
    r["solvent_conditions"]["composition"] = F("MeCN", ["MeCN"], evidence)
    r["solvent_conditions"]["concentration"] = F("0.1 M", Q(0.1, "M"), evidence)
    r["solvent_conditions"]["concentration_basis"] = F("MeCN (0.1 M)", "nominal reaction concentration", evidence)


table4_yields = dict(zip(list("abcdefghijklmnopqrst"), [80,77,67,89,90,49,90,67,89,56,76,52,79,67,80,66,67,95,90,64]))
for letter, y in table4_yields.items():
    if letter <= "k":
        aa_mid = "m_1a"
        allene_mid = ensure(f"m_terminal_allene_4{letter}", f"terminal allene for 4{letter}", 6,
                            note=f"Identity retained from the Table 4 source drawing for product 4{letter}; no complete SMILES asserted.")
    else:
        aa_mid = ensure(f"m_aminoacid_4{letter}", f"amino acid substrate for 4{letter}", 6,
                        note=f"Identity retained from the Table 4 amino-acid source drawing for product 4{letter}.")
        allene_mid = "m_terminal_allene_oh"
    pmid = ensure(f"m_4{letter}", f"4{letter}", 6, stereo_kind="product_linear")
    r = base_reaction(f"table4_4{letter}", 6, "specific_experiment", ["substrate_scope"], f"Table 4 product 4{letter}",
                      "decarboxylative addition of amino acids to terminal allenes", "linear decarboxylative hydroaminoalkylation")
    terminal_scope_participants(r, ev(6), aa_mid, allene_mid)
    r["participants"].append(part("product", pmid, "product", ev(6), product_role="main product"))
    reaction_step(r, ev(6), [p["participant_id"] for p in r["participants"]], "irradiation under Pd/rac-BINAP/photoredox conditions",
                  ("16 h", Q(16, "h"), False), ("rt", "room temperature", False))
    r["special_conditions"]["light_source"] = F("blue LED", "blue LED", ev(6))
    outcome(r, "yield", "yield", [pmid], f"{y}% yield", Q(y, "%"), ev(6), "isolated yield", yield_type="isolated yield")
    if letter == "k":
        outcome(r, "lb", "selectivity", [pmid], "77:23 L/B", "77:23", ev(6), None, scope="linear/branched ratio")
        outcome(r, "ez", "selectivity", [pmid], ">95:5 E/Z", ">95:5", ev(6), None, scope="E/Z ratio")
    else:
        outcome(r, "lb", "selectivity", [pmid], ">95:5 L/B (unless otherwise specified)", ">95:5", ev(6), None, scope="linear/branched ratio")
        outcome(r, "ez", "selectivity", [pmid], "≫95:5 E/Z (unless otherwise specified)", "≫95:5", ev(6), None, scope="E/Z ratio")
    r["classification"] = class_partial(ev(6), ["C-C", "C-H"])
    REACTIONS.append(r)


# Scheme 2 derivatizations.
for label, kind in [("5a", "specified_drawing"), ("5b", "specified_drawing"), ("5c", "product_branched"), ("5d", "specified_drawing")]:
    ensure(f"m_{label}", label, 6, stereo_kind=kind)

def derivatization(rid, start_mid, product_mid, steps, y, eratio, labels):
    r = base_reaction(rid, 6, "specific_experiment", ["derivatization"], f"Scheme 2 {product_mid.replace('m_','')}",
                      "synthetic application of allylic adduct", "postcatalysis derivatization")
    evidence = ev(6)
    r["participants"] = [part("start", start_mid, "reactant", evidence), part("product", product_mid, "product", evidence, product_role="main product")]
    pids = ["start", "product"]
    for i, spec in enumerate(steps, 1):
        for pid, mid, role, kwargs in spec["participants"]:
            r["participants"].append(part(pid, mid, role, evidence, **kwargs))
            pids.append(pid)
        s = st.new_step(f"reaction_{i}", evidence)
        s["participant_ids"] = [x[0] for x in spec["participants"]] + (["start"] if i == 1 else []) + (["product"] if i == len(steps) else [])
        s["operation"] = F(spec["operation"], spec["operation"], evidence)
        if spec.get("temp"): s["temperature"] = F(spec["temp"], spec["temp"], evidence)
        if spec.get("pressure"): s["pressure"] = F(spec["pressure"][0], spec["pressure"][1], evidence)
        r["steps"].append(s)
    outcome(r, "yield", "yield", [product_mid], f"{y}% yield", Q(y, "%"), evidence, "yield basis not stated in Scheme 2", yield_type="reported yield")
    outcome(r, "er", "selectivity", [product_mid], f"{eratio} er", eratio, evidence, None, scope="enantiomeric ratio")
    r["classification"] = class_partial(evidence, labels)
    r["review"]["unresolved"].append("Scheme 2 does not report step times, workup/purification, or the yield/er measurement methods.")
    REACTIONS.append(r)

derivatization("scheme2_5a", "m_3al", "m_5a", [{"operation":"TFA/DCM deprotection at RT", "temp":"room temperature", "participants":[("tfa","m_tfa","reagent",{}),("dcm","m_dcm","solvent",{})]}], 85, "99:1", ["O-H"])
derivatization("scheme2_5b", "m_3ag", "m_5b", [
    {"operation":"hydrogenation with Pd/C in MeOH", "temp":"room temperature", "pressure":("H2 (1 atm)", "H2 at 1 atm"), "participants":[("pdc","m_pdc","catalyst",{"ref":"m_3ag","molpct":("10%",Q(10,"mol%"))}),("h2","m_h2","reagent",{}),("meoh","m_meoh","solvent",{})]},
    {"operation":"TFA/DCM deprotection", "temp":"room temperature", "participants":[("tfa","m_tfa","reagent",{}),("dcm","m_dcm","solvent",{})]}], 85, "91.5:8.5", ["C-H","O-H"])
derivatization("scheme2_5c", "m_3aa", "m_5c", [{"operation":"hydroboration then oxidation", "participants":[("bbn","m_9bbn","reagent",{}),("h2o2","m_h2o2","reagent",{})]}], 63, "96:4", ["C-O","C-H"])
derivatization("scheme2_5d", "m_3aa", "m_5d", [
    {"operation":"Rh-catalyzed hydroformylation", "temp":{"number":65,"unit":"°C"}, "pressure":("CO/H2 (1/1; 6 atm)", "CO/H2 1:1 at 6 atm total pressure"), "participants":[("rh","m_rh","catalyst",{"ref":"m_3aa","molpct":("2 mol %",Q(2,"mol%"))}),("biphephos","m_biphephos","ligand",{"ref":"m_3aa","molpct":("4 mol %",Q(4,"mol%"))}),("co","m_co","reagent",{}),("h2","m_h2","reagent",{}),("thf","m_thf","solvent",{})]},
    {"operation":"hydrogenation in MeOH/EtOAc 1/1", "pressure":("H2 (1 atm)", "H2 at 1 atm"), "participants":[("pdc","m_pdc","catalyst",{"ref":"m_3aa","molpct":("10 mol %",Q(10,"mol%"))}),("h2b","m_h2","reagent",{}),("meoh","m_meoh","solvent",{}),("etoac","m_etoac","solvent",{})]}], 60, "96:4", ["C-C","C-N","C-H"])


# Mechanistic experiments reported in pages 4–5 prose. Sparse records deliberately preserve missing SI details.
def mech(rid, page, title, reactant_mids, agents, outcome_kind, raw_outcome, labels=None, product_mid=None):
    r = base_reaction(rid, page, "specific_experiment", ["mechanistic_study"], f"page {page} mechanistic-study prose", title, "mechanistic probe")
    evidence = ev(page)
    for i, mid in enumerate(reactant_mids): r["participants"].append(part(f"r{i+1}", mid, "reactant", evidence))
    for i, mid in enumerate(agents): r["participants"].append(part(f"a{i+1}", mid, "reagent", evidence))
    reaction_step(r, evidence, [p["participant_id"] for p in r["participants"]], title)
    targets = []
    if product_mid:
        r["participants"].append(part("product", product_mid, "product", evidence, product_role="observed product"))
        targets = [product_mid]
    elif outcome_kind == "selectivity":
        r["participants"].append(part("product", "m_probe_product", "product", evidence, product_role="probe product"))
        targets = ["m_probe_product"]
    outcome(r, "observation", outcome_kind, targets, raw_outcome, raw_outcome, evidence, "method/detail not stated in main-text prose")
    r["classification"] = class_partial(evidence, labels) if labels else class_na(evidence, "Analytical/control observation does not establish a net synthetic bond change.")
    r["review"]["unresolved"].append("Amounts, exact conditions, and measurement details are deferred to unavailable Supporting Information.")
    REACTIONS.append(r)

mech("mech_tempo_inhibition",4,"TEMPO radical-scavenger control",["m_1a","m_2a"],["m_tempo"],"no_reaction","reaction completely inhibited")
ensure("m_giese_product", "α-alkylated Giese-type product", 4, note="Product mentioned but structure/label not printed in the covered prose.")
mech("mech_ethyl_acrylate_trap",4,"ethyl acrylate radical-trap experiment",["m_1a"],["m_ethyl_acrylate"],"other","α-alkylated product formed through a Giese-type reaction",["C-C"],"m_giese_product")
ensure("m_probe_product", "hydroaminoalkylation probe product", 4, stereo_kind="product_branched")
mech("mech_L1a_probe",4,"stereochemical probe using L-1a",["m_L1a","m_2a"],[],"selectivity","same regio- and stereoselectivities as D-1a",["C-C","C-H"])
mech("mech_D1a_probe",4,"stereochemical probe using D-1a",["m_D1a","m_2a"],[],"selectivity","same regio- and stereoselectivities as L-1a",["C-C","C-H"])
ensure("m_1a_dimer", "dimer compound from 1a", 4, note="Dimer observed in mechanistic experiment; structure not provided on covered pages.")
mech("mech_cyclohexylallene_dimer",4,"coupling of 1a with cyclohexyl allene under linear-selective conditions",["m_1a","m_cyclohexyl_allene"],[],"other","dimer compound from 1a was produced",["C-C"],"m_1a_dimer")
mech("mech_stern_volmer_1a",5,"Stern–Volmer study of 1a with PC1",["m_1a"],["m_pc1"],"other","1a was not able to quench photoexcited Ir(III)")
mech("mech_stern_volmer_cs2co3",5,"Stern–Volmer study of basified 1a with PC1",["m_1a"],["m_cs2co3","m_pc1"],"other","effective reductive quenching observed")
mech("mech_quantum_yield",5,"quantum-yield measurement",["m_1a","m_2a"],["m_pc1"],"other","low quantum yield; radical-chain propagation not involved")
mech("mech_light_on_off",5,"light on/off experiment",["m_1a","m_2a"],["m_pc1"],"other","irradiation required throughout the reaction")
ensure("m_5aa", "5aa", 5, name="ester intermediate 5aa", note="Ester intermediate observed by 1H NMR; full structure not encoded.")
mech("mech_nmr_ester_5aa",5,"1H NMR monitoring during the first hour",["m_1a","m_2a"],[],"other","fast ester 5aa formation observed within the first 1 h",["C-O"],"m_5aa")

# Page-5 prose mechanism: retained as proposed, not an experimentally proven sequence.
r = base_reaction("proposed_mechanism_text",5,"proposed_mechanism",[],"page 5 prose describing Scheme 3",
                  "proposed mechanism for decarboxylative hydroaminoalkylation",None)
r["classification"] = class_unresolved(ev(5), "Proposed elementary intermediates are not a net concrete experiment and Scheme 3 itself is outside the bounded page images.")
r["review"]["unresolved"].append("Only the page-5 prose description is in the bounded review; the Scheme 3 drawing is on page 7 and is not claimed as covered.")
REACTIONS.append(r)


# Assemble and save the formal bundle.
bundle = st.read(RUN / "submitted.json")
bundle["molecules"] = list(MOLECULES.values())
bundle["reactions"] = REACTIONS
bundle["source_context"] = [
    {"context_id": f"{DOC}:article_metadata", "document_id": DOC, "kind": "article_metadata",
     "title": "Dual Palladium/Photoredox-Catalyzed Enantioselective and Regioselective Decarboxylative Hydroaminoalkylation of Allenes",
     "authors": ["Jun Zheng", "Ali Nikbakht", "Bernhard Breit"], "doi": "10.1021/acscatal.1c00153", "evidence_ids": ev(1)},
    {"context_id": f"{DOC}:page1_graphical_abstract", "document_id": DOC, "kind": "generic_summary",
     "note": "Graphical abstract summarizes branched and linear reaction families; it is not an additional experiment.", "evidence_ids": ev(1)},
    {"context_id": f"{DOC}:figure1_bioactives", "document_id": DOC, "kind": "background_examples",
     "note": "Figure 1 lists bioactive 1,2-amino-alcohol examples; these are not experiments performed in this article.", "evidence_ids": ev(2)},
    {"context_id": f"{DOC}:si_status", "document_id": DOC, "kind": "source_limit",
     "note": "Supporting Information is mentioned in the article but was not provided or explicitly associated in this workspace.", "evidence_ids": ev(1)+ev(4)+ev(5)},
]
st.write(RUN / "submitted.json", bundle)

raw = {
    "document_id": DOC, "pages": [1,2,3,4,5,6], "source": "assistant visual transcription plus extracted-text cross-check",
    "table1_rows": opt_rows, "table2_rows": table2, "table3_rows": table3,
    "table3_special": [{"product":"3al","yield":83,"ligand":"ent-L9","dr":["7.3:1",">19:1","12:1"],"B/L":"91:9"},
                       {"product":"3al'","yield":60,"ligand":"L9","dr":["7.1:1","2.2:1","2:1"],"B/L":"91:9"}],
    "table4_yields": table4_yields,
    "deduplication": {"3aa": "Table 3 3aa merged as additional evidence into Table 1 entry 2"},
    "semantic_accuracy_verified": False,
}
st.write(RUN / "raw_model/pages_0001_0006_extraction.json", raw)

page_records = {i: [] for i in range(1,7)}
for rr in REACTIONS:
    for eid in rr["evidence_ids"]:
        for page in range(1,7):
            if f":page{page}:" in eid and rr["reaction_id"] not in page_records[page]:
                page_records[page].append(rr["reaction_id"])
coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": [
    {"page":1,"state":"attempted","regions":[
        {"region":"title, abstract, and graphical abstract","finding":"Article metadata and a generic two-branch reaction summary; no independent concrete experiment.","evidence_ids":ev(1)},
        {"region":"body text","finding":"Background and cited prior work; not converted into current-paper experiments.","evidence_ids":ev(1)}],
     "record_ids":page_records[1],"unresolved":["Graphical structures are generic and were not treated as concrete molecule instances."]},
    {"page":2,"state":"attempted","regions":[
        {"region":"Figure 1","finding":"Representative bioactive compounds supplied as background, not synthetic experiments in this paper.","evidence_ids":ev(2)},
        {"region":"Scheme 1A","finding":"Conceptual retrosynthetic analysis.","evidence_ids":ev(2)},
        {"region":"Scheme 1B","finding":"Cited prior-work generic alkoxyallene allylations, classified as literature_example.","evidence_ids":ev(2)},
        {"region":"Scheme 1C","finding":"Generic scheme for this work, not a separate concrete experiment.","evidence_ids":ev(2)}],
     "record_ids":page_records[2],"unresolved":["Generic R groups prevent complete molecule identities and atom mapping."]},
    {"page":3,"state":"attempted","regions":[
        {"region":"optimization prose","finding":"One initial (S)-BINAP experiment and discussion of optimization results.","evidence_ids":ev(3)},
        {"region":"Table 1 rows 1–8","finding":"Five productive optimization rows and three no-reaction controls, each recorded separately.","evidence_ids":ev(3)},
        {"region":"Table 1 footnote","finding":"General charges, 16 h/rt/blue LED, NMR yield basis, isolated-yield parentheses, crude-1H-NMR dr, and chiral-HPLC er.","evidence_ids":ev(3)}],
     "record_ids":page_records[3],"unresolved":["No SI procedure, workup, purification, atmosphere, LED wavelength or power was supplied."]},
    {"page":4,"state":"attempted","regions":[
        {"region":"Table 2 header and footnotes","finding":"Standard branched conditions and measurement methods; PC1 and derivatization exceptions retained by row.","evidence_ids":ev(4)},
        {"region":"Table 2 products 3ba–3pa","finding":"Fifteen distinct amino-acid scope experiments transcribed row by row.","evidence_ids":ev(4)},
        {"region":"lower-column prose","finding":"Table 4 summary, Scheme 2 summary, and five real mechanistic probes (TEMPO, trap, L/D probes, dimer) identified.","evidence_ids":ev(4)}],
     "record_ids":page_records[4],"unresolved":["Most complex drawn structures were retained by source label rather than guessed SMILES; detailed mechanistic conditions are in unavailable SI."]},
    {"page":5,"state":"attempted","regions":[
        {"region":"Table 3 header and products","finding":"Fourteen displayed results inspected; 3aa is duplicate evidence for Table 1 entry 2, while 13 other condition/product rows are distinct.","evidence_ids":ev(5)},
        {"region":"mechanistic prose","finding":"Two Stern–Volmer observations, low quantum yield, light on/off experiment, and first-hour 1H-NMR ester observation recorded separately.","evidence_ids":ev(5)},
        {"region":"proposed mechanism prose","finding":"Author-proposed catalytic sequence classified as proposed_mechanism, not experimental proof.","evidence_ids":ev(5)}],
     "record_ids":page_records[5],"unresolved":["Scheme 3 drawing is on page 7 and is outside this bounded visual review; 3al/3al' absolute configurations are not assigned on page 5."]},
    {"page":6,"state":"attempted","regions":[
        {"region":"Table 4 products 4a–4t","finding":"Twenty terminal-allene scope rows transcribed separately with yields and default or exceptional L/B and E/Z ratios.","evidence_ids":ev(6)},
        {"region":"Table 4 footnote","finding":"0.2 mmol scale, 1.2 equiv allene, Pd/rac-BINAP/PC1 loadings, MeCN 0.1 M, 16 h, rt, blue LED, isolated-yield basis and default selectivities.","evidence_ids":ev(6)},
        {"region":"Scheme 2 four arrows","finding":"Four derivatizations (5a–5d), including the two-step 5b and 5d sequences, recorded without converting overall yield into a last-step yield.","evidence_ids":ev(6)}],
     "record_ids":page_records[6],"unresolved":["Scheme 2 omits times, workup/purification, and yield/er measurement methods; several substrate/product graphs remain source-drawing identities."]}
]}
st.write(RUN / "coverage/pages_0001_0006.json", coverage)

mentions = st.source_mentions(st.read(RUN / "dataset.json"), ["1a","2a","3aa","3al","4a","5a","5aa"], 180)
st.write(RUN / "raw_model/source_mentions_pages_0001_0006.json", mentions)

print(json.dumps({"molecules":len(bundle["molecules"]),"reactions":len(bundle["reactions"]),
                  "coverage_pages":[p["page"] for p in coverage["pages"]]}, ensure_ascii=False, indent=2))
