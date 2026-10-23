import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_d12c5bb9ecbba832"
E = {n: [f"{DOC}:page{n}:text", f"{DOC}:page{n}:image"] for n in range(1, 7)}


def fact(raw=None, normalized=None, evidence=(), method="original_extraction",
         status=None, derivation=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def derived(raw, normalized, evidence, why):
    return fact(raw, normalized, evidence, method="structure_derivation", derivation=why)


def inferred(raw, normalized, evidence, why):
    return fact(raw, normalized, evidence, method="inference", derivation=why)


def missing(status="unreported"):
    return fact(status=status)


def source_molecule(mid, label, name, evidence, smiles=None, stereo=None, note=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = fact(label, label, evidence)
    m["names"]["original"] = fact(name, name, evidence)
    m["names"]["english"] = fact(name, name, evidence)
    m["names"]["abbreviations"] = fact(label, [label], evidence) if label else missing()
    m["names"]["aliases"] = missing()
    m["cas_numbers"] = missing()
    stereo = stereo or {}
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            "candidate connectivity interpreted from the source name/drawing",
            smiles,
            evidence,
            "Candidate connectivity was transcribed from the source name and structure drawing; parsing does not verify the source identity.",
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing()
    for key in ("atom_centers", "double_bonds", "axial", "mixture", "assignment"):
        if key in stereo:
            raw, normalized = stereo[key]
            m["stereochemistry"][key] = fact(raw, normalized, evidence)
        else:
            m["stereochemistry"][key] = missing()
    for key, default in (("chirality_status", "unresolved"),
                         ("geometry_status", "unresolved"),
                         ("representation_status", "unresolved")):
        val = stereo.get(key, default)
        if isinstance(val, tuple):
            m["stereochemistry"][key] = fact(val[0], val[1], evidence)
        else:
            m["stereochemistry"][key] = inferred(
                "stereochemical review state",
                val,
                evidence,
                "Status records what the supplied source specifies and whether the candidate SMILES preserves it; unresolved does not mean achiral or racemic.",
            )
    if smiles and stereo.get("encode_graph"):
        diag = st.stereo_inspect(smiles)
        m["stereochemistry"]["isomeric_smiles"] = derived(
            "tool-normalized candidate isomeric SMILES",
            diag["canonical_isomeric_smiles"],
            evidence,
            "Generated from the candidate SMILES with standard stereo diagnostics; only the source-specified alkene geometry is encoded.",
        )
        m["stereochemistry"]["identity_key"] = derived(
            "stereo graph key for the candidate representation",
            diag["stereo_graph_key"],
            evidence,
            "SHA-256 key of the canonical candidate isomeric graph; this is not a sample identity or source verification.",
        )
    else:
        m["stereochemistry"]["isomeric_smiles"] = missing()
        m["stereochemistry"]["identity_key"] = missing()
    m["identity_note"] = fact(
        note or "Machine-created molecule entry; source name/label retained and candidate structure remains subject to expert review.",
        note or "Machine-created molecule entry; source name/label retained and candidate structure remains subject to expert review.",
        evidence,
    )
    m["external_properties"] = []
    return m


def named_unresolved(mid, name, evidence, label=None):
    return source_molecule(
        mid, label or name, name, evidence, smiles=None,
        note="Named participant is source-linked, but no concrete structure was asserted in this extraction because the paper did not supply a structure drawing for this participant.",
    )


def generic_molecule(mid, label, name, evidence, fixed, variables):
    m = source_molecule(mid, label, name, evidence, smiles=None,
                        note="Generic source drawing contains variable substituents; it is retained as a template and not expanded to a concrete molecule.")
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": fact(fixed, {"description": fixed}, evidence),
        "variables": [],
        "instance_mapping": [],
        "reason": fact("generic structure in source figure", "Source figure contains variable or generic substituents.", evidence),
    }
    for var in variables:
        m["generic_template"]["variables"].append({
            "label": var,
            "attachment_points": fact(f"attachment of {var} shown in source drawing", f"source-drawn attachment for {var}", evidence),
            "definition": fact(f"{var} is not uniquely defined in this figure", None, evidence, status="unreadable"),
            "allowed_values": missing(),
            "constraints": missing(),
        })
    return m


def reaction(rid, record_type, evidence, location, categories, name, transformation, unresolved):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = fact(", ".join(categories) if categories else "none", categories, evidence)
    r["location"] = fact(location, {"description": location}, evidence)
    r["reaction_name"] = fact(name, name, evidence) if name else missing()
    r["transformation_type"] = fact(transformation, transformation, evidence) if transformation else missing()
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved)}
    r["completeness"] = "partial"
    return r


def participant(r, pid, mid, role, evidence, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    p["product_role"] = fact(product_role, product_role, evidence) if product_role else missing()
    r["participants"].append(p)
    return p


def charge(p, kind, raw, normalized, evidence, reference=None):
    p["charges"][kind] = fact(raw, normalized, evidence)
    if kind in ("equivalents", "mol_percent"):
        p["equivalence_reference"] = fact(reference, reference, evidence)


def scale(r, mid, raw, normalized, evidence):
    r["reaction_scale"] = fact(raw, {"reference_molecule_id": mid, "amount": normalized}, evidence)


def solvents(r, raw, mids, evidence, total=None):
    r["solvent_conditions"]["composition"] = fact(raw, mids, evidence)
    r["solvent_conditions"]["total_volume"] = fact(raw, total, evidence) if total else missing()


def step(r, sid, stage, pids, evidence, operation, temperature=None, time=None,
         temperature_fact=None, atmosphere=True):
    s = st.new_step(sid, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = fact(operation, operation, evidence)
    s["temperature"] = temperature_fact if temperature_fact else (
        fact(temperature[0], temperature[1], evidence) if temperature else missing())
    s["time"] = fact(time[0], time[1], evidence) if time else missing()
    s["atmosphere"] = fact(
        "All reactions were carried out ... under an atmosphere of argon unless otherwise specified.",
        "argon", E[3], method="inherited_general_condition",
    ) if atmosphere and stage == "reaction" else missing()
    s["pressure"] = missing()
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, evidence, raw, normalized, product_role=None,
            yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["product_role"] = fact(product_role, product_role, evidence) if product_role else missing()
    o["value"] = fact(raw, normalized, evidence)
    o["yield_type"] = fact(yield_type, yield_type, evidence) if yield_type else missing()
    o["measurement_basis"] = fact(basis, basis, evidence) if basis else missing()
    o["scope"] = fact(scope, scope, evidence) if scope else missing()
    r["outcomes"].append(o)
    return o


def procedure(r, raw, evidence, addition=None, reference=None):
    r["procedure"]["original_text"] = fact(raw, raw, evidence)
    if addition:
        r["procedure"]["addition_order"] = fact(addition, addition, evidence)
    if reference:
        r["procedure"]["general_procedure_reference"] = fact(reference, reference, evidence)


def note(r, raw, evidence, normalized=None, method="original_extraction", derivation=None):
    r["notes"].append(fact(raw, normalized if normalized is not None else raw,
                           evidence, method=method, derivation=derivation))


def set_classification(r, evidence, changes, limitation):
    normalized_changes = []
    labels = {}
    for i, item in enumerate(changes, 1):
        a, b, before, after, kind, scope_text = item
        cid = f"b{i}"
        ch = {"change_id": cid,
              "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
              "before_order": before, "after_order": after, "change_type": kind,
              "scope": scope_text}
        normalized_changes.append(ch)
        if kind in ("formed", "hydrogen_gain", "order_increased"):
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], dict(lab, supporting_change_ids=[]))["supporting_change_ids"].append(cid)
    r["classification"]["mapped_reaction_smiles"] = missing()
    r["classification"]["bond_changes"] = [
        derived(f"graph change inferred from source drawing: {c['scope']}", c, evidence,
                "Bond change was manually inferred from the source drawing without asserting complete atom mapping.")
        for c in normalized_changes
    ]
    r["classification"]["labels"] = derived(
        "bond-formation taxonomy labels inferred from the listed partial graph changes",
        list(labels.values()), evidence,
        "Labels include all and only formed/promoted bonds in the partial change list; atom mapping and omitted coproducts remain unresolved.",
    )
    r["classification"]["coverage"] = inferred(
        "partial graph comparison", "partial", evidence,
        "Only source-visible connectivity changes were encoded; no complete atom mapping or balance is claimed.",
    )
    r["classification"]["limitations"] = [limitation, "No three-model votes were generated."]
    r["classification"]["model_votes"] = []


dataset = json.loads((RUN / "dataset.json").read_text())
initialized_bundle = json.loads((RUN / "submitted.json").read_text())
assert initialized_bundle["schema_version"] == "3.1"
bundle = st.empty(copy.deepcopy(dataset["documents"]))

# Source-defined principal compounds.
bundle["molecules"].extend([
    source_molecule("m_epoxyalkyne_h", "trans-2-ethynyl-3-isopropyloxirane",
                    "trans-2-ethynyl-3-isopropyloxirane", E[3], "CC(C)C1OC1C#C",
                    {"assignment": ("trans-2-ethynyl-3-isopropyloxirane", {"ring_relationship": "trans", "absolute_configuration": "unassigned"}),
                     "chirality_status": ("trans; absolute configuration not reported", "partially_specified"),
                     "geometry_status": ("oxirane trans relationship", "specified"),
                     "representation_status": ("candidate SMILES omits relative ring stereochemistry", "not_encoded")}),
    source_molecule("m7", "7", "trans-2-(deuterioethynyl)-3-isopropyloxirane", E[3] + E[4],
                    "CC(C)C1OC1C#C[2H]",
                    {"assignment": ("Trans-2-(deuterioethynyl)-3-isopropyloxirane", {"ring_relationship": "trans", "absolute_configuration": "unassigned"}),
                     "chirality_status": ("trans; absolute configuration not reported", "partially_specified"),
                     "geometry_status": ("oxirane trans relationship", "specified"),
                     "representation_status": ("candidate SMILES omits relative ring stereochemistry", "not_encoded")}),
    source_molecule("m8", "8", "(3R*,5R*)-6-deuterio-2-methylhexa-4,5-dien-3-ol", E[2] + E[4],
                    "CC(C)C(O)C=C=[CH][2H]",
                    {"atom_centers": ("(3R*,5R*)", [{"position": 3, "descriptor": "R*", "assignment": "relative"}]),
                     "axial": ("(3R*,5R*) allene assignment", [{"allene_axis": "C4=C5=C6", "descriptor": "5R*", "assignment": "relative"}]),
                     "assignment": ("apparently one predominant stereoisomer; confirmed retrospectively", "relative configuration assigned; absolute configuration not reported"),
                     "chirality_status": ("(3R*,5R*)", "specified"),
                     "geometry_status": ("allene axial stereochemistry is not alkene E/Z geometry", "not_applicable"),
                     "representation_status": ("candidate SMILES does not encode the source relative center/allene axis", "not_encoded")}),
    source_molecule("m8_hydroxycarbamate", "8-hydroxycarbamate",
                    "(3R*,5R*)-6-deuterio-2-methylhexa-4,5-dien-3-yl hydroxycarbamate", E[4],
                    "CC(C)C(OC(=O)NO)C=C=[CH][2H]",
                    {"atom_centers": ("(3R*,5R*)", [{"position": 3, "descriptor": "R*", "assignment": "relative"}]),
                     "axial": ("(3R*,5R*)", [{"allene_axis": "C4=C5=C6", "descriptor": "5R*", "assignment": "relative"}]),
                     "chirality_status": ("(3R*,5R*)", "specified"), "geometry_status": ("allene axial stereochemistry is not alkene E/Z geometry", "not_applicable"),
                     "representation_status": ("relative point/axial stereochemistry not encoded", "not_encoded")}),
    source_molecule("m9", "9", "(3R*,5R*)-6-deuterio-2-methylhexa-4,5-dien-3-yl p-toluenesulfonyloxycarbamate", E[2] + E[4],
                    "CC(C)C(OC(=O)NOS(=O)(=O)c1ccc(C)cc1)C=C=[CH][2H]",
                    {"atom_centers": ("(3R*,5R*)", [{"position": 3, "descriptor": "R*", "assignment": "relative"}]),
                     "axial": ("(3R*,5R*)", [{"allene_axis": "C4=C5=C6", "descriptor": "5R*", "assignment": "relative"}]),
                     "chirality_status": ("(3R*,5R*)", "specified"), "geometry_status": ("allene axial stereochemistry is not alkene E/Z geometry", "not_applicable"),
                     "representation_status": ("relative point/axial stereochemistry not encoded", "not_encoded")}),
    source_molecule("m10", "10", "(3R*,5R*)-6-deuterio-2-methylhexa-4,5-dien-3-yl carbamate", E[2],
                    "CC(C)C(OC(=O)N)C=C=[CH][2H]",
                    {"atom_centers": ("Scheme 3 carries the stereodefined allene from 8 to 10", [{"position": 3, "descriptor": "R*", "assignment": "relative"}]),
                     "axial": ("Scheme 3 stereodefined allene", [{"allene_axis": "C4=C5=C6", "descriptor": "5R*", "assignment": "relative"}]),
                     "chirality_status": ("relative configuration inherited in Scheme 3", "specified"), "geometry_status": ("allene axial stereochemistry is not alkene E/Z geometry", "not_applicable"),
                     "representation_status": ("relative point/axial stereochemistry not encoded", "not_encoded")}),
    source_molecule("m4", "4", "(4R*,5S*,Z)-4-isopropyl-6-(methylene-d)-3-oxa-1-azabicyclo[3.1.0]hexane-2-one", E[2] + E[3] + E[4],
                    "CC(C)C1OC(=O)N2C1/C2=C/[2H]",
                    {"atom_centers": ("(4R*,5S*)", [{"position": 4, "descriptor": "R*"}, {"position": 5, "descriptor": "S*"}]),
                     "double_bonds": ("(Z)-6-(methylene-d)", [{"bond": "C6=methylene-d", "configuration": "Z"}]),
                     "assignment": ("NOE correlations and 3JHH = 9.5 Hz support the assigned stereochemistry", "relative center and alkene assignments supported by NMR/NOE; absolute configuration not reported"),
                     "chirality_status": ("(4R*,5S*)", "specified"), "geometry_status": ("Z", "specified"),
                     "representation_status": ("Z geometry encoded; relative ring stereochemistry not encoded", "partial"), "encode_graph": True}),
    source_molecule("m11", "11", "(4S*,5R*)-5-isopropyl-4-[(E)-propen-2-yl-1-d]oxazolidin-2-one", E[2] + E[3] + E[4] + E[5],
                    "CC(C)C1OC(=O)NC1/C(C)=C/[2H]",
                    {"atom_centers": ("(4S*,5R*)", [{"position": 4, "descriptor": "S*"}, {"position": 5, "descriptor": "R*"}]),
                     "double_bonds": ("[(E)-propen-2-yl-1-d]", [{"bond": "propenyl C=C", "configuration": "E"}]),
                     "assignment": ("methyl group cis to deuterium; diagnostic NOE correlations", "relative stereochemistry supported by NMR/NOE; absolute configuration not reported"),
                     "chirality_status": ("(4S*,5R*)", "specified"), "geometry_status": ("E", "specified"),
                     "representation_status": ("E geometry encoded; relative ring stereochemistry not encoded", "partial"), "encode_graph": True}),
    source_molecule("m12", "12", "(4S*,5R*)-4-[(E)-buta-1,3-dien-2-yl-1-d]-5-isopropyloxazolidin-2-one", E[2] + E[3] + E[5],
                    "CC(C)C1OC(=O)NC1/C(C=C)=C/[2H]",
                    {"atom_centers": ("(4S*,5R*)", [{"position": 4, "descriptor": "S*"}, {"position": 5, "descriptor": "R*"}]),
                     "double_bonds": ("[(E)-buta-1,3-dien-2-yl-1-d]", [{"bond": "diene C1=C2", "configuration": "E"}, {"bond": "terminal vinyl C=C", "configuration": "unspecified/not stereogenic"}]),
                     "assignment": ("new C-C bond cis to deuterium; diagnostic NOE correlations", "relative stereochemistry supported by NMR/NOE; absolute configuration not reported"),
                     "chirality_status": ("(4S*,5R*)", "specified"), "geometry_status": ("E", "specified"),
                     "representation_status": ("E geometry encoded; relative ring stereochemistry not encoded", "partial"), "encode_graph": True}),
    source_molecule("m13", "13", "(R*)-2-methyl-1-[(S*)-3-(methyl-d)-2H-azirin-2-yl]propyl pent-4-enoate", E[2] + E[3] + E[5],
                    "C=CCCC(=O)OC(C(C)C)C1N=C1C[2H]",
                    {"atom_centers": ("(R*) ... [(S*) ...]", [{"site": "2-methyl-1-substituted propyl center", "descriptor": "R*"}, {"site": "azirine C2", "descriptor": "S*"}]),
                     "assignment": ("(R*)/(S*) relative descriptors", "relative configuration stated; absolute configuration not reported"),
                     "chirality_status": ("(R*)/(S*)", "specified"), "geometry_status": ("terminal alkene", "not_applicable"),
                     "representation_status": ("relative stereochemistry not encoded", "not_encoded")}),
])

# Generic structures drawn only to frame prior work or proposed mechanisms.
for args in [
    ("m_scheme1_1", "1", "generic allene carbamate 1", E[1], "Scheme 1 allene carbamate framework", ["R′"]),
    ("m_scheme1_2", "2", "generic bicyclic methyleneaziridine 2", E[1], "Scheme 1 fused oxazolidinone methyleneaziridine framework", ["R′"]),
    ("m_scheme1_3", "3", "generic S_NV ring-opened oxazolidinone 3", E[1], "Scheme 1 ring-opened oxazolidinone framework", ["R′"]),
    ("m5", "5", "generic cis-Nu/D product 5", E[2], "Scheme 2 retentive cis-Nu/D oxazolidinone product", ["Nu"]),
    ("m6", "6", "generic trans-Nu/D product 6", E[2], "Scheme 2 invertive trans-Nu/D oxazolidinone product", ["Nu"]),
    ("m14", "14", "generic Cu(III) intermediate 14", E[3], "Scheme 5 formal C-N insertion Cu intermediate", ["R", "M"]),
    ("m15", "15", "generic organocopper ligand-coupling product 15", E[3], "Scheme 5 ligand-coupled N-Cu product", ["R"]),
]:
    bundle["molecules"].append(generic_molecule(*args))

# Named reagents, solvents, and workup/purification agents are kept source-linked but structurally unresolved.
named = {
    "m_k2co3": ("K2CO3", E[2] + E[3] + E[4]), "m_me_cn": ("acetonitrile", E[2] + E[3] + E[4]),
    "m_d2o": ("D2O", E[2] + E[3]), "m_petrol": ("petrol (30–40 °C petroleum ether fraction)", E[3] + E[4] + E[5]),
    "m_mgso4": ("MgSO4", E[3] + E[4]), "m_dibal": ("DIBAL", E[2] + E[4]),
    "m_dcm": ("dichloromethane", E[2] + E[4]), "m_water": ("water", E[4]),
    "m_rochelle": ("saturated aqueous Rochelle's salt", E[4]), "m_brine": ("brine", E[4] + E[5]),
    "m_cdi": ("N,N-carbonyldiimidazole", E[2] + E[4]), "m_imidazole": ("imidazole", E[4]),
    "m_nh2oh_hcl": ("NH2OH·HCl", E[2] + E[4]), "m_hcl": ("hydrochloric acid", E[4]),
    "m_ether": ("diethyl ether", E[2] + E[3] + E[4] + E[5]), "m_na2so4": ("Na2SO4", E[4] + E[5]),
    "m_tscl": ("p-TsCl", E[2] + E[4]), "m_et3n": ("triethylamine", E[2] + E[4]),
    "m_rh2oac4": ("Rh2(OAc)4", E[2] + E[4]), "m_meli": ("methyllithium", E[3] + E[4]),
    "m_cui": ("CuI", E[1] + E[3] + E[4] + E[5]), "m_thf": ("THF", E[3] + E[4] + E[5]),
    "m_nh4cl": ("saturated aqueous NH4Cl", E[4] + E[5]), "m_vinylmgbr": ("vinylmagnesium bromide", E[2] + E[3] + E[5]),
    "m_trichloro_isocyanate": ("trichloroacetyl isocyanate", E[2]), "m_meoh": ("methanol", E[2]),
    "m_phioac2": ("PhI(OAc)2", E[2]), "m_mgo": ("MgO", E[2]),
    "m_rhii_generic": ("Rh(II) catalyst (Rh2(OAc)4, Rh2(esp)2, or Rh2(TPA)4)", E[2]),
    "m_prior_organo": ("R′MgX/CuI or Ph2CuLi", E[1]), "m_generic_nu": ("generic nucleophile Nu−", E[2]),
    "m_generic_cuprate": ("generic R2Cu− M+", E[3]),
}
for mid, (name, ev) in named.items():
    bundle["molecules"].append(named_unresolved(mid, name, ev))

# Scheme 1: cited prior work, not experiments newly reported in this paper.
r = reaction("scheme1_prior_aziridination", "literature_example", E[1], "page 1, Scheme 1, first arrow",
             [], "intramolecular allene aziridination (prior work, Ref. 4)", "allene aziridination",
             ["Conditions and yield are delegated to Ref. 4 and are not reported here."])
participant(r, "s1_1", "m_scheme1_1", "reactant", E[1]); participant(r, "s1_2", "m_scheme1_2", "product", E[1], product_role="generic product")
step(r, "scheme_arrow", "reaction", ["s1_1", "s1_2"], E[1], "Ref. 4", atmosphere=False)
procedure(r, "Ref. 4", E[1], reference="Ref. 4")
set_classification(r, E[1], [("C", "N", 0, 1, "formed", "first aziridine C-N bond"), ("C", "N", 0, 1, "formed", "second aziridine C-N bond")],
                   "Generic prior-work structures and absent conditions prevent complete atom mapping.")
bundle["reactions"].append(r)

r = reaction("scheme1_prior_snv", "literature_example", E[1], "page 1, Scheme 1, second arrow",
             [], "S_NV ring-opening (prior work)", "nucleophilic vinylic substitution",
             ["The figure summarizes prior work; individual R′MgX examples and yields are not enumerated in this paper."])
participant(r, "s1_2", "m_scheme1_2", "reactant", E[1]); participant(r, "s1_reagent", "m_prior_organo", "reagent", E[1]); participant(r, "s1_3", "m_scheme1_3", "product", E[1], product_role="generic product")
step(r, "scheme_arrow", "reaction", ["s1_2", "s1_reagent", "s1_3"], E[1], "R′MgX/[CuI] or Ph2CuLi", atmosphere=False)
procedure(r, "R′MgX, [CuI] or Ph2CuLi; [R = i-Pr]", E[1])
set_classification(r, E[1], [("C", "N", 1, 0, "broken", "exo C-N cleavage"), ("C", "C", 0, 1, "formed", "new carbon-carbon bond to R′"), ("N", "H", 0, 1, "hydrogen_gain", "protonated oxazolidinone nitrogen")],
                   "Prior-work generic substituents prevent complete mapping.")
bundle["reactions"].append(r)

# Scheme 2 proposed alternatives, explicitly not observed experiments.
for rid, product, mode, relationship in [
    ("scheme2_retention_path", "m5", "S_NVπ", "cis-Nu/D"),
    ("scheme2_inversion_path", "m6", "S_NVσ", "trans-Nu/D"),
]:
    r = reaction(rid, "proposed_mechanism", E[2], f"page 2, Scheme 2, proposed {mode} path",
                 [], f"proposed {mode} pathway", "proposed nucleophilic vinylic substitution mechanism",
                 ["Generic Nu− and transition-state drawings are mechanistic proposals, not isolated experiments."])
    participant(r, "substrate4", "m4", "reactant", E[2]); participant(r, "nu", "m_generic_nu", "reagent", E[2]); participant(r, "hypothetical_product", product, "product", E[2], product_role="hypothetical mechanistic product")
    step(r, "proposed_path", "reaction", ["substrate4", "nu", "hypothetical_product"], E[2], f"proposed {mode} path leading to {relationship}", atmosphere=False)
    procedure(r, f"Simplified {mode} reaction mode of methyleneaziridine 4 with generic nucleophile Nu−", E[2])
    set_classification(r, E[2], [("C", "N", 1, 0, "broken", "proposed exo C-N cleavage"), ("C", "C", 0, 1, "formed", "proposed Nu carbon attachment"), ("N", "H", 0, 1, "hydrogen_gain", "proposed protonation shown after H+")],
                       "Hypothetical generic mechanism; bonds are schematic and not a balanced experimental mapping.")
    bundle["reactions"].append(r)

# Experimental deuteration to 7.
r = reaction("exp_deuteration_to_7", "specific_experiment", E[3] + E[4], "pages 3–4, experimental section for compound 7",
             ["mechanistic_study"], "H/D exchange to trans-2-(deuterioethynyl)-3-isopropyloxirane", "isotopic exchange",
             ["Temperature and deuterium-incorporation measurement method are not reported.", "The candidate structures retain the trans assignment textually but do not encode an absolute enantiomer."])
p0 = participant(r, "epoxyalkyne", "m_epoxyalkyne_h", "reactant", E[3]); charge(p0, "mass", "2.03 g", {"number": 2.03, "unit": "g"}, E[3]); charge(p0, "amount", "18.4 mmol", {"number": 18.4, "unit": "mmol"}, E[3])
pk = participant(r, "k2co3", "m_k2co3", "reagent", E[3]); charge(pk, "mass", "3.76 g", {"number": 3.76, "unit": "g"}, E[3]); charge(pk, "amount", "27.2 mmol", {"number": 27.2, "unit": "mmol"}, E[3])
ps = participant(r, "mecn", "m_me_cn", "solvent", E[3]); charge(ps, "volume", "42 mL", {"number": 42, "unit": "mL"}, E[3])
pd = participant(r, "d2o", "m_d2o", "reagent", E[3]); charge(pd, "volume", "20 mL", {"number": 20, "unit": "mL"}, E[3])
participant(r, "product7", "m7", "product", E[3] + E[4], product_role="main product")
pp = participant(r, "petrol_workup", "m_petrol", "workup_agent", E[3], stage="workup"); charge(pp, "volume", "5 × 100 mL", {"number": 100, "unit": "mL", "count": 5}, E[3])
participant(r, "mgso4", "m_mgso4", "workup_agent", E[4], stage="workup")
scale(r, "m_epoxyalkyne_h", "18.4 mmol trans-2-ethynyl-3-isopropyloxirane", {"number": 18.4, "unit": "mmol"}, E[3])
solvents(r, "acetonitrile (42 mL)", ["m_me_cn"], E[3], {"number": 42, "unit": "mL"})
step(r, "pre_stir", "reaction", ["epoxyalkyne", "k2co3", "mecn"], E[3], "stir before adding D2O", time=("30 min", {"number": 30, "unit": "min"}))
step(r, "exchange", "reaction", ["epoxyalkyne", "k2co3", "d2o", "mecn"], E[3], "add D2O and continue stirring", time=("5 h", {"number": 5, "unit": "h"}))
step(r, "extract", "workup", ["petrol_workup"], E[3], "extract into petrol (5 × 100 mL)", atmosphere=False)
step(r, "dry_and_concentrate", "workup", ["mgso4"], E[3] + E[4], "dry combined extracts with MgSO4 and remove solvent in vacuo; product is volatile", atmosphere=False)
procedure(r, "Trans-2-ethynyl-3-isopropyloxirane was added to K2CO3 in acetonitrile; after 30 min D2O was added and stirring continued for 5 h.", E[3], addition="epoxyalkyne to K2CO3/acetonitrile; D2O after 30 min")
outcome(r, "yield7", "yield", ["m7"], E[4], "73%", {"number": 73, "unit": "%"}, "main product", "reported yield", "isolated after extraction/drying/solvent removal")
outcome(r, "amount7", "product_amount", ["m7"], E[4], "1.50 g", {"number": 1.50, "unit": "g"}, "main product", basis="isolated material")
outcome(r, "d_incorp7", "selectivity", ["m7"], E[4], "94% deuterium incorporation", {"number": 94, "unit": "%"}, scope="deuterium incorporation at the ethynyl terminus")
set_classification(r, E[3] + E[4], [("C", "H", 1, 0, "hydrogen_loss", "terminal protium removed"), ("C", "H", 0, 1, "formed", "terminal C-D bond formed; D classified as hydrogen isotope")],
                   "Isotopic exchange is encoded without complete atom/isotope mapping or a stated proton/deuteron-transfer mechanism.")
bundle["reactions"].append(r)

# 7 -> 8 coordinated reductive ring opening.
r = reaction("exp_7_to_8", "specific_experiment", E[2] + E[4], "page 2 Scheme 3(i) and page 4 experimental section for 8",
             ["mechanistic_study"], "coordinated reductive ring-opening of an alkynyl epoxide", "alkynyl epoxide to stereodefined allene alcohol",
             ["Addition duration, quench-water volume, reaction concentration, and absolute configuration are unreported.", "Hydride/proton atom sources and full atom mapping remain unresolved."])
p7 = participant(r, "substrate7", "m7", "reactant", E[4]); charge(p7, "mass", "1.50 g", {"number": 1.50, "unit": "g"}, E[4]); charge(p7, "amount", "13.5 mmol", {"number": 13.5, "unit": "mmol"}, E[4])
pdibal = participant(r, "dibal", "m_dibal", "reagent", E[2] + E[4]); charge(pdibal, "volume", "33.7 mL of a 1.0 M solution in hexane", {"number": 33.7, "unit": "mL", "solution_concentration": {"number": 1.0, "unit": "M"}, "carrier": "hexane"}, E[4]); charge(pdibal, "amount", "33.7 mmol", {"number": 33.7, "unit": "mmol"}, E[4])
pdcm1 = participant(r, "dcm_substrate", "m_dcm", "solvent", E[4]); charge(pdcm1, "volume", "100 mL", {"number": 100, "unit": "mL"}, E[4])
pdcm2 = participant(r, "dcm_dibal", "m_dcm", "solvent", E[4]); charge(pdcm2, "volume", "100 mL", {"number": 100, "unit": "mL"}, E[4])
participant(r, "product8", "m8", "product", E[2] + E[4], product_role="main product")
participant(r, "water_quench", "m_water", "workup_agent", E[4], stage="workup")
pr = participant(r, "rochelle", "m_rochelle", "workup_agent", E[4], stage="workup"); charge(pr, "volume", "200 mL", {"number": 200, "unit": "mL"}, E[4])
pex = participant(r, "dcm_extract", "m_dcm", "workup_agent", E[4], stage="workup"); charge(pex, "volume", "5 × 50 mL", {"number": 50, "unit": "mL", "count": 5}, E[4])
pb = participant(r, "brine_wash", "m_brine", "workup_agent", E[4], stage="workup"); charge(pb, "volume", "50 mL", {"number": 50, "unit": "mL"}, E[4])
participant(r, "mgso4", "m_mgso4", "workup_agent", E[4], stage="workup")
participant(r, "petrol_purification", "m_petrol", "purification_agent", E[4], stage="purification"); participant(r, "ether_purification", "m_ether", "purification_agent", E[4], stage="purification")
scale(r, "m7", "13.5 mmol epoxyalkyne 7", {"number": 13.5, "unit": "mmol"}, E[4])
solvents(r, "dichloromethane: 100 mL substrate solution added to 100 mL DIBAL solution", ["m_dcm"], E[4])
step(r, "dropwise_addition", "reaction", ["substrate7", "dibal", "dcm_substrate", "dcm_dibal"], E[4], "add solution of 7 dropwise to DIBAL solution", temperature=("0 °C", {"number": 0, "unit": "°C"}))
step(r, "hold", "reaction", ["substrate7", "dibal"], E[2] + E[4], "stir reaction mixture", temperature=("0 °C", {"number": 0, "unit": "°C"}), time=("1 h", {"number": 1, "unit": "h"}))
step(r, "quench", "workup", ["water_quench", "rochelle"], E[4], "quench carefully with water; add saturated Rochelle's salt dropwise and stir overnight", time=("overnight", "overnight"), atmosphere=False)
step(r, "extract_wash_dry", "workup", ["dcm_extract", "brine_wash", "mgso4"], E[4], "extract with dichloromethane, wash with brine, dry with MgSO4, and remove solvent", atmosphere=False)
step(r, "flash", "purification", ["petrol_purification", "ether_purification"], E[4], "flash chromatography (petrol/ether, 8:1)", atmosphere=False)
procedure(r, "A solution of epoxyalkyne 7 in dichloromethane was added dropwise to DIBAL in dichloromethane at 0 °C; the mixture was stirred for 1 h and worked up with water and Rochelle's salt.", E[4], addition="7 solution into DIBAL solution at 0 °C")
outcome(r, "yield8", "yield", ["m8"], E[2] + E[4], "53%", {"number": 53, "unit": "%"}, "main product", "reported yield", "isolated after flash chromatography")
outcome(r, "amount8", "product_amount", ["m8"], E[4], "810 mg", {"number": 810, "unit": "mg"}, "main product", basis="isolated material")
note(r, "allene 8 was isolated apparently as one predominant stereoisomer, depicted as that expected, and confirmed retrospectively from the NMR data for methyleneaziridine 4", E[2])
set_classification(r, E[2] + E[4], [("C", "O", 1, 0, "broken", "epoxide C-O bond cleavage"), ("C", "C", 1, 2, "order_increased", "allene C=C bond formed from the epoxide-to-alkyne connecting bond"), ("C", "C", 3, 2, "order_decreased", "alkyne bond reduced to allene double bond"), ("C", "H", 0, 1, "hydrogen_gain", "hydride-derived C-H bond"), ("O", "H", 0, 1, "hydrogen_gain", "alcohol O-H after workup")],
                   "Connectivity changes are inferred from Scheme 3; hydride/proton provenance and complete atom mapping are not established.")
bundle["reactions"].append(r)

# 8 -> hydroxycarbamate.
r = reaction("exp_8_to_hydroxycarbamate", "specific_experiment", E[2] + E[4], "page 4, first stage under compound 9",
             ["derivatization"], "formation of hydroxycarbamate from allene alcohol 8", "carbamoylation",
             ["The paper reports approximate TLC completion times but no reaction temperature other than RT and no concentration basis.", "Yield measurement is isolated after flash chromatography; the basis of deuterium/stereochemical integrity is not separately quantified."])
p8 = participant(r, "alcohol8", "m8", "reactant", E[4]); charge(p8, "mass", "500 mg", {"number": 500, "unit": "mg"}, E[4]); charge(p8, "amount", "4.42 mmol", {"number": 4.42, "unit": "mmol"}, E[4])
pcdi = participant(r, "cdi", "m_cdi", "reagent", E[4]); charge(pcdi, "mass", "1.43 g", {"number": 1.43, "unit": "g"}, E[4]); charge(pcdi, "amount", "8.82 mmol", {"number": 8.82, "unit": "mmol"}, E[4])
pmec = participant(r, "mecn", "m_me_cn", "solvent", E[4]); charge(pmec, "volume", "25 mL", {"number": 25, "unit": "mL"}, E[4])
pim = participant(r, "imidazole", "m_imidazole", "reagent", E[4]); charge(pim, "mass", "1.20 g", {"number": 1.20, "unit": "g"}, E[4]); charge(pim, "amount", "17.6 mmol", {"number": 17.6, "unit": "mmol"}, E[4])
pnh = participant(r, "nh2oh_hcl", "m_nh2oh_hcl", "reagent", E[4]); charge(pnh, "mass", "1.54 g", {"number": 1.54, "unit": "g"}, E[4]); charge(pnh, "amount", "22.2 mmol", {"number": 22.2, "unit": "mmol"}, E[4])
participant(r, "hydroxycarbamate", "m8_hydroxycarbamate", "product", E[4], product_role="main product")
phcl = participant(r, "hcl_partition", "m_hcl", "workup_agent", E[4], stage="workup"); charge(phcl, "volume", "12 mL, 1.0 M", {"number": 12, "unit": "mL", "solution_concentration": {"number": 1.0, "unit": "M"}}, E[4])
pe = participant(r, "ether_partition", "m_ether", "workup_agent", E[4], stage="workup"); charge(pe, "volume", "10 mL", {"number": 10, "unit": "mL"}, E[4])
pex = participant(r, "ether_extractions", "m_ether", "workup_agent", E[4], stage="workup"); charge(pex, "volume", "3 × 50 mL", {"number": 50, "unit": "mL", "count": 3}, E[4])
pbr = participant(r, "brine", "m_brine", "workup_agent", E[4], stage="workup"); charge(pbr, "volume", "25 mL", {"number": 25, "unit": "mL"}, E[4])
participant(r, "na2so4", "m_na2so4", "workup_agent", E[4], stage="workup"); participant(r, "petrol_purification", "m_petrol", "purification_agent", E[4], stage="purification"); participant(r, "ether_purification", "m_ether", "purification_agent", E[4], stage="purification")
scale(r, "m8", "4.42 mmol alcohol 8", {"number": 4.42, "unit": "mmol"}, E[4]); solvents(r, "acetonitrile (25 mL)", ["m_me_cn"], E[4], {"number": 25, "unit": "mL"})
step(r, "cdi_stage", "reaction", ["alcohol8", "cdi", "mecn"], E[2] + E[4], "add CDI and stir until complete by TLC", temperature=("RT", "r.t."), time=("~24 h", {"number": 24, "unit": "h", "qualifier": "approx", "raw_text": "~24 h"}))
step(r, "hydroxylamine_stage", "reaction", ["imidazole", "nh2oh_hcl"], E[2] + E[4], "add imidazole and NH2OH·HCl and continue stirring until complete by TLC", temperature=("RT", "r.t."), time=("~24 h", {"number": 24, "unit": "h", "qualifier": "approx", "raw_text": "~24 h"}))
step(r, "partition_extract", "workup", ["hcl_partition", "ether_partition", "ether_extractions", "brine", "na2so4"], E[4], "filter, remove solvent, partition between HCl and ether, extract, wash with brine, dry, and remove solvent", atmosphere=False)
step(r, "flash", "purification", ["petrol_purification", "ether_purification"], E[4], "flash chromatography (petrol/ether, 3:1 to 1:1)", atmosphere=False)
procedure(r, "CDI was added to alcohol 8 in acetonitrile and stirred at RT until complete by TLC (~24 h); imidazole and NH2OH·HCl were then added and stirring continued (~24 h).", E[4], addition="CDI first; imidazole and NH2OH·HCl after TLC completion")
outcome(r, "yield_hydroxycarbamate", "yield", ["m8_hydroxycarbamate"], E[4], "60%", {"number": 60, "unit": "%"}, "main product", "isolated yield", "isolated after flash chromatography")
outcome(r, "brsm_yield_hydroxycarbamate", "yield", ["m8_hydroxycarbamate"], E[4], "73% based on recovered 8", {"number": 73, "unit": "%"}, "main product", "yield based on recovered starting material", "based on recovered alcohol 8")
outcome(r, "amount_hydroxycarbamate", "product_amount", ["m8_hydroxycarbamate"], E[4], "454 mg", {"number": 454, "unit": "mg"}, "main product", basis="isolated material")
outcome(r, "recovered_8", "recovery", ["m8"], E[4], "90 mg, 18%", {"number": 18, "unit": "%", "mass": {"number": 90, "unit": "mg"}}, yield_type="recovered starting material", basis="recovered after flash chromatography")
set_classification(r, E[2] + E[4], [("C", "O", 0, 1, "formed", "alcohol oxygen bonded to carbamoyl carbon")],
                   "Carbamoyl atom sources and leaving groups are not completely mapped.")
bundle["reactions"].append(r)

# Hydroxycarbamate -> 9.
r = reaction("exp_hydroxycarbamate_to_9", "specific_experiment", E[2] + E[4], "page 2 Scheme 3(ii)(b) and page 4 second stage under compound 9",
             ["derivatization"], "tosylation of hydroxycarbamate", "O-sulfonylation",
             ["The scheme gives 0 °C to RT over 18 h, while the prose specifies addition at 0 °C and only states that stirring continued for 18 h; exact ramp/hold subdivision is unreported."])
ph = participant(r, "hydroxycarbamate", "m8_hydroxycarbamate", "reactant", E[4]); charge(ph, "mass", "554 mg", {"number": 554, "unit": "mg"}, E[4]); charge(ph, "amount", "3.22 mmol", {"number": 3.22, "unit": "mmol"}, E[4])
pts = participant(r, "tscl", "m_tscl", "reagent", E[4]); charge(pts, "mass", "619 mg", {"number": 619, "unit": "mg"}, E[4]); charge(pts, "amount", "3.25 mmol", {"number": 3.25, "unit": "mmol"}, E[4])
pet = participant(r, "ether_solvent", "m_ether", "solvent", E[4]); charge(pet, "volume", "30 mL", {"number": 30, "unit": "mL"}, E[4])
pte = participant(r, "et3n", "m_et3n", "reagent", E[4]); charge(pte, "volume", "0.45 mL", {"number": 0.45, "unit": "mL"}, E[4]); charge(pte, "amount", "3.23 mmol", {"number": 3.23, "unit": "mmol"}, E[4])
participant(r, "product9", "m9", "product", E[4], product_role="main product")
ped = participant(r, "ether_dilution", "m_ether", "workup_agent", E[4], stage="workup"); charge(ped, "volume", "30 mL", {"number": 30, "unit": "mL"}, E[4])
pbr = participant(r, "brine", "m_brine", "workup_agent", E[4], stage="workup"); charge(pbr, "volume", "2 × 20 mL", {"number": 20, "unit": "mL", "count": 2}, E[4])
participant(r, "na2so4", "m_na2so4", "workup_agent", E[4], stage="workup"); participant(r, "petrol_purification", "m_petrol", "purification_agent", E[4], stage="purification"); participant(r, "ether_purification", "m_ether", "purification_agent", E[4], stage="purification")
scale(r, "m8_hydroxycarbamate", "3.22 mmol hydroxycarbamate", {"number": 3.22, "unit": "mmol"}, E[4]); solvents(r, "dry ether (30 mL)", ["m_ether"], E[4], {"number": 30, "unit": "mL"})
step(r, "et3n_addition", "reaction", ["hydroxycarbamate", "tscl", "ether_solvent", "et3n"], E[4], "add triethylamine dropwise to hydroxycarbamate and p-TsCl", temperature=("0 °C", {"number": 0, "unit": "°C"}))
step(r, "warm_and_stir", "reaction", ["hydroxycarbamate", "tscl", "et3n"], E[2] + E[4], "continue stirring while warming from 0 °C to RT", temperature_fact=fact("0 °C to RT", None, E[2] + E[4], status="pending_review"), time=("18 h", {"number": 18, "unit": "h"}))
step(r, "wash_dry", "workup", ["ether_dilution", "brine", "na2so4"], E[4], "dilute with ether, wash with brine, dry with Na2SO4, and remove solvent", atmosphere=False)
step(r, "flash", "purification", ["petrol_purification", "ether_purification"], E[4], "flash chromatography (petrol/ether, 5:1 to pure ether)", atmosphere=False)
procedure(r, "Recrystallised p-TsCl was added to the hydroxycarbamate in dry ether at 0 °C; triethylamine was added dropwise and stirring continued for 18 h.", E[4], addition="p-TsCl with hydroxycarbamate; triethylamine dropwise at 0 °C")
outcome(r, "yield9", "yield", ["m9"], E[4], "85%", {"number": 85, "unit": "%"}, "main product", "isolated yield", "isolated after flash chromatography")
outcome(r, "amount9", "product_amount", ["m9"], E[4], "893 mg", {"number": 893, "unit": "mg"}, "main product", basis="isolated material")
set_classification(r, E[2] + E[4], [("O", "S", 0, 1, "formed", "hydroxycarbamate N-O oxygen bonded to tosyl sulfur")],
                   "Leaving groups and acid/base coproducts are not mapped.")
bundle["reactions"].append(r)

# Scheme-only 8 -> 10.
r = reaction("scheme3_8_to_10", "specific_experiment", E[2], "page 2, Scheme 3(iii)", ["derivatization"],
             "formation of carbamate 10", "carbamoylation",
             ["This experiment is only summarized in Scheme 3; all charges, scale, workup, purification, and yield measurement basis are unreported."])
participant(r, "alcohol8", "m8", "reactant", E[2]); participant(r, "trichloroacetyl_isocyanate", "m_trichloro_isocyanate", "reagent", E[2]); participant(r, "dcm", "m_dcm", "solvent", E[2]); participant(r, "k2co3", "m_k2co3", "reagent", E[2]); participant(r, "meoh", "m_meoh", "solvent", E[2]); participant(r, "product10", "m10", "product", E[2], product_role="main product")
solvents(r, "CH2Cl2, then MeOH", ["m_dcm", "m_meoh"], E[2])
step(r, "isocyanate_stage", "reaction", ["alcohol8", "trichloroacetyl_isocyanate", "dcm"], E[2], "treat with trichloroacetyl isocyanate in CH2Cl2", temperature=("0 °C", {"number": 0, "unit": "°C"}), time=("4 h", {"number": 4, "unit": "h"}))
step(r, "methanolysis_stage", "reaction", ["k2co3", "meoh"], E[2], "then treat with K2CO3 in MeOH", temperature=("RT", "r.t."), time=("4 h", {"number": 4, "unit": "h"}))
procedure(r, "Cl3CCO·NCO, CH2Cl2, 0 °C, 4 h then K2CO3, MeOH, RT, 4 h (10, 97%)", E[2])
outcome(r, "yield10", "yield", ["m10"], E[2], "97%", {"number": 97, "unit": "%"}, "main product", "reported scheme yield")
set_classification(r, E[2], [("C", "O", 0, 1, "formed", "alcohol oxygen bonded to carbamoyl carbon")],
                   "Scheme-only record lacks charges, workup, purification, and full atom mapping.")
bundle["reactions"].append(r)

# 9 -> 4 tethered allene aziridination.
r = reaction("exp_9_to_4", "specific_experiment", E[2] + E[4], "page 2 Scheme 3(iv) and page 4 experimental section for 4",
             ["mechanistic_study"], "diastereoselective tethered allene aziridination", "intramolecular allene aziridination",
             ["Absolute configuration is not assigned; the candidate SMILES encodes Z geometry but not the relative ring centers.", "Concentration and catalyst-preparation details are unreported."])
p9 = participant(r, "carbamate9", "m9", "reactant", E[4]); charge(p9, "mass", "100 mg", {"number": 100, "unit": "mg"}, E[4]); charge(p9, "amount", "0.306 mmol", {"number": 0.306, "unit": "mmol"}, E[4])
prh = participant(r, "rh2oac4", "m_rh2oac4", "catalyst", E[2] + E[4]); charge(prh, "mass", "6.8 mg", {"number": 6.8, "unit": "mg"}, E[4]); charge(prh, "amount", "15.4 μmol", {"number": 15.4, "unit": "umol"}, E[4]); charge(prh, "mol_percent", "5 mol%", {"number": 5, "unit": "mol%"}, E[2], reference="m9")
pk = participant(r, "k2co3", "m_k2co3", "reagent", E[2] + E[4]); charge(pk, "mass", "127 mg", {"number": 127, "unit": "mg"}, E[4]); charge(pk, "amount", "0.919 mmol", {"number": 0.919, "unit": "mmol"}, E[4])
pm = participant(r, "mecn", "m_me_cn", "solvent", E[2] + E[4]); charge(pm, "volume", "3 mL", {"number": 3, "unit": "mL"}, E[4])
participant(r, "product4", "m4", "product", E[2] + E[4], product_role="main product")
participant(r, "mecn_dilution", "m_me_cn", "workup_agent", E[4], stage="workup"); participant(r, "petrol_purification", "m_petrol", "purification_agent", E[4], stage="purification"); participant(r, "ether_purification", "m_ether", "purification_agent", E[4], stage="purification")
scale(r, "m9", "0.306 mmol carbamate 9", {"number": 0.306, "unit": "mmol"}, E[4]); solvents(r, "dry acetonitrile (3 mL)", ["m_me_cn"], E[4], {"number": 3, "unit": "mL"})
step(r, "aziridination", "reaction", ["carbamate9", "rh2oac4", "k2co3", "mecn"], E[2] + E[4], "add Rh2(OAc)4 and K2CO3 and stir vigorously", temperature=("25 °C", {"number": 25, "unit": "°C"}), time=("90 min", {"number": 90, "unit": "min"}))
step(r, "filter_concentrate", "workup", ["mecn_dilution"], E[4], "dilute with acetonitrile, filter, and concentrate in vacuo at RT", temperature=("RT", "r.t."), atmosphere=False)
step(r, "flash", "purification", ["petrol_purification", "ether_purification"], E[4], "flash chromatography (petrol/ether, 7:1)", atmosphere=False)
procedure(r, "Rh2(OAc)4 and K2CO3 were added to carbamate 9 in dry acetonitrile at 25 °C; the mixture was stirred vigorously for 90 min.", E[4], addition="Rh2(OAc)4 and K2CO3 added to the solution of 9")
outcome(r, "yield4", "yield", ["m4"], E[2] + E[4], "25%", {"number": 25, "unit": "%"}, "main product", "isolated yield", "isolated after flash chromatography")
outcome(r, "amount4", "product_amount", ["m4"], E[4], "12 mg", {"number": 12, "unit": "mg"}, "main product", basis="isolated material")
note(r, "The stereochemistry in methyleneaziridine 4 was confirmed by comparisons with non-deuterated 2 and the NOE correlations shown in Fig. 1.", E[2])
set_classification(r, E[2] + E[4], [("C", "N", 0, 1, "formed", "first bond of the tethered aziridine ring"), ("C", "N", 0, 1, "formed", "second bond of the tethered aziridine ring"), ("C", "C", 2, 1, "order_decreased", "allene double bond consumed in aziridination")],
                   "Two C-N bond formations are inferred from the source-drawn intramolecular aziridination; complete atom mapping is not asserted.")
bundle["reactions"].append(r)

# Scheme 3(v): catalyst family cannot be decomposed into individual experiments from the paper.
r = reaction("scheme3_10_to_4_rh_screen", "generic_scheme", E[2], "page 2 Scheme 3(v) and adjacent discussion",
             ["optimization"], "Rh(II)-catalysed allene aziridination from carbamate 10", "intramolecular allene aziridination",
             ["A range of Rh(II) catalysts and yields is reported without catalyst-to-yield pairing; individual experiments cannot be reconstructed.", "All charges, scale, workup, purification, and yield measurement basis are unreported."])
participant(r, "carbamate10", "m10", "reactant", E[2]); prh = participant(r, "rhii", "m_rhii_generic", "catalyst", E[2]); charge(prh, "mol_percent", "5 mol%", {"number": 5, "unit": "mol%"}, E[2], reference="m10")
participant(r, "phioac2", "m_phioac2", "reagent", E[2]); participant(r, "mgo", "m_mgo", "additive", E[2]); participant(r, "dcm", "m_dcm", "solvent", E[2]); participant(r, "product4", "m4", "product", E[2], product_role="reported product across unresolved catalyst screen")
solvents(r, "CH2Cl2", ["m_dcm"], E[2]); step(r, "screen_family", "reaction", ["carbamate10", "rhii", "phioac2", "mgo", "dcm"], E[2], "Rh(II) catalyst/PhI(OAc)2/MgO in CH2Cl2", temperature=("reflux", "reflux"), time=("48 h", {"number": 48, "unit": "h"}))
procedure(r, "[from 10] Rh(II) catalyst (5 mol%; see text), PhI(OAc)2, MgO, CH2Cl2, reflux, 48 h (4–17%)", E[2])
outcome(r, "yield4_range", "yield", ["m4"], E[2], "4–17%", {"minimum": 4, "maximum": 17, "unit": "%", "raw_text": "4–17%"}, "reported product across unresolved catalyst screen", "range across catalyst variants", scope="unresolved set of Rh2(OAc)4, Rh2(esp)2, and Rh2(TPA)4 experiments")
set_classification(r, E[2], [("C", "N", 0, 1, "formed", "first bond of the tethered aziridine ring"), ("C", "N", 0, 1, "formed", "second bond of the tethered aziridine ring"), ("C", "C", 2, 1, "order_decreased", "allene double bond consumed in aziridination")],
                   "This generic record represents an unresolved catalyst screen, not one concrete experiment.")
bundle["reactions"].append(r)

# 4 -> 11. Preserve scheme/prose temperature conflict.
r = reaction("exp_4_to_11", "specific_experiment", E[2] + E[3] + E[4] + E[5], "page 3 Scheme 4(i), pages 4–5 experimental section for 11",
             ["mechanistic_study"], "S_NV ring-opening with lithium dimethylcuprate", "stereoretentive nucleophilic vinylic substitution",
             ["Temperature program conflicts: Scheme 4 gives −20 °C to RT over 30 min; prose places warming to RT after the quench.", "Relative ring stereochemistry is not encoded in the candidate SMILES; only the E alkene geometry is encoded."])
pmeli = participant(r, "meli", "m_meli", "reagent", E[4]); charge(pmeli, "volume", "0.16 mL of a 1.6 M solution in ether", {"number": 0.16, "unit": "mL", "solution_concentration": {"number": 1.6, "unit": "M"}, "carrier": "diethyl ether"}, E[4]); charge(pmeli, "amount", "0.256 mmol", {"number": 0.256, "unit": "mmol"}, E[4])
pcui = participant(r, "cui", "m_cui", "reagent", E[4]); charge(pcui, "mass", "25 mg", {"number": 25, "unit": "mg"}, E[4]); charge(pcui, "amount", "0.131 mmol", {"number": 0.131, "unit": "mmol"}, E[4])
pthf1 = participant(r, "thf_cui", "m_thf", "solvent", E[4]); charge(pthf1, "volume", "1 mL", {"number": 1, "unit": "mL"}, E[4])
p4 = participant(r, "substrate4", "m4", "reactant", E[4]); charge(p4, "mass", "20 mg", {"number": 20, "unit": "mg"}, E[4]); charge(p4, "amount", "0.130 mmol", {"number": 0.130, "unit": "mmol"}, E[4])
pthf2 = participant(r, "thf_substrate", "m_thf", "solvent", E[4]); charge(pthf2, "volume", "1 mL", {"number": 1, "unit": "mL"}, E[4])
participant(r, "product11", "m11", "product", E[2] + E[3] + E[4] + E[5], product_role="main product")
pnh = participant(r, "nh4cl", "m_nh4cl", "workup_agent", E[4], stage="workup"); charge(pnh, "volume", "3 mL", {"number": 3, "unit": "mL"}, E[4])
pex = participant(r, "ether_extract", "m_ether", "workup_agent", E[4], stage="workup"); charge(pex, "volume", "3 × 20 mL", {"number": 20, "unit": "mL", "count": 3}, E[4])
pbr = participant(r, "brine", "m_brine", "workup_agent", E[4], stage="workup"); charge(pbr, "volume", "5 mL", {"number": 5, "unit": "mL"}, E[4])
participant(r, "na2so4", "m_na2so4", "workup_agent", E[4], stage="workup"); participant(r, "ether_recrystallization", "m_ether", "purification_agent", E[4], stage="purification")
scale(r, "m4", "0.130 mmol methyleneaziridine 4", {"number": 0.130, "unit": "mmol"}, E[4]); solvents(r, "THF: 1 mL CuI suspension and 1 mL substrate solution; methyllithium supplied in ether", ["m_thf", "m_ether"], E[4])
step(r, "cuprate_preparation", "reaction", ["meli", "cui", "thf_cui"], E[4], "add methyllithium dropwise to CuI suspension and wait", temperature=("−20 °C", {"number": -20, "unit": "°C"}), time=("15 min", {"number": 15, "unit": "min"}))
temp_conflict = fact(["Scheme 4: −20 °C → RT, 30 min", "prose: stirred 30 min, quenched, then allowed to warm to RT"], None, E[3] + E[4], status="conflict", alternatives=[fact("−20 °C → RT", "−20 °C to r.t.", E[3]), fact("warming to RT after quench", "warming after quench", E[4])])
step(r, "substitution", "reaction", ["substrate4", "meli", "cui", "thf_substrate"], E[3] + E[4], "add solution of 4 and stir for 30 min", temperature_fact=temp_conflict, time=("30 min", {"number": 30, "unit": "min"}))
step(r, "quench_and_warm", "workup", ["nh4cl"], E[4], "quench with saturated aqueous NH4Cl and then allow to warm to RT", temperature=("warm to RT", "r.t."), atmosphere=False)
step(r, "extract_wash_dry", "workup", ["ether_extract", "brine", "na2so4"], E[4], "extract with ether, wash with brine, dry with Na2SO4, and concentrate", atmosphere=False)
step(r, "recrystallize", "purification", ["ether_recrystallization"], E[4], "recrystallize from ether", atmosphere=False)
procedure(r, "Methyllithium was added dropwise to CuI in THF at −20 °C; after 15 min a solution of 4 in THF was added and the mixture was stirred for 30 min before NH4Cl quench.", E[4], addition="MeLi to CuI; after 15 min add 4 solution")
outcome(r, "yield11", "yield", ["m11"], E[2] + E[3] + E[4], "77%", {"number": 77, "unit": "%"}, "main product", "isolated yield", "isolated after recrystallization")
outcome(r, "amount11", "product_amount", ["m11"], E[4], "17 mg", {"number": 17, "unit": "mg"}, "main product", basis="isolated material")
outcome(r, "stereo11", "selectivity", ["m11"], E[2] + E[3], "methyl group cis to the deuterium atom", "cis-to-D; stereoretentive", "main product", basis="NMR experiments including diagnostic NOE", scope="relative disposition of newly formed C-C bond/methyl substituent and deuterium")
note(r, "In 11 no NOE correlation was observed between the vinyl methyl protons and =CHD.", E[2] + E[3])
set_classification(r, E[2] + E[3] + E[4], [("C", "N", 1, 0, "broken", "exo vinylic C-N bond cleavage"), ("C", "C", 0, 1, "formed", "new vinylic carbon-methyl bond"), ("N", "H", 0, 1, "hydrogen_gain", "oxazolidinone N-H after workup")],
                   "Methyl origin is source-linked, but complete organocopper/proton atom mapping and coproducts are unresolved.")
bundle["reactions"].append(r)

# 4 -> 12 + 13, one experimental operation with two isolated products.
r = reaction("exp_4_to_12_13", "specific_experiment", E[2] + E[3] + E[5], "page 3 Scheme 4(ii) and page 5 experimental section for 12/13",
             ["mechanistic_study"], "copper-catalysed Grignard S_NV ring-opening with competing azirine formation", "stereoretentive nucleophilic vinylic substitution with side reaction",
             ["Temperature program conflicts: Scheme 4 gives −50 °C to 0 °C over 1 h; prose gives −50 °C to −20 °C over 1 h before quench.", "The exact atom mapping and measurement basis for the proposed carbonyl-addition/1,4-vinylation route to 13 are unresolved.", "Relative ring stereochemistry is not encoded in the candidate SMILES; only the E alkene geometry of 12 is encoded."])
pvmg = participant(r, "vinylmgbr", "m_vinylmgbr", "reagent", E[3] + E[5]); charge(pvmg, "volume", "0.25 mL of a 1.0 M solution in THF", {"number": 0.25, "unit": "mL", "solution_concentration": {"number": 1.0, "unit": "M"}, "carrier": "THF"}, E[5]); charge(pvmg, "amount", "0.25 mmol", {"number": 0.25, "unit": "mmol"}, E[5]); charge(pvmg, "equivalents", "2.0 eq", {"number": 2.0, "unit": "equiv"}, E[3], reference="m4")
pcui = participant(r, "cui", "m_cui", "catalyst", E[3] + E[5]); charge(pcui, "mass", "1.2 mg", {"number": 1.2, "unit": "mg"}, E[5]); charge(pcui, "amount", "6.30 μmol", {"number": 6.30, "unit": "umol"}, E[5]); charge(pcui, "mol_percent", "5 mol%", {"number": 5, "unit": "mol%"}, E[3], reference="m4")
pthf1 = participant(r, "thf_cui", "m_thf", "solvent", E[5]); charge(pthf1, "volume", "1 mL", {"number": 1, "unit": "mL"}, E[5])
p4 = participant(r, "substrate4", "m4", "reactant", E[5]); charge(p4, "mass", "19 mg", {"number": 19, "unit": "mg"}, E[5]); charge(p4, "amount", "0.123 mmol", {"number": 0.123, "unit": "mmol"}, E[5])
pthf2 = participant(r, "thf_substrate", "m_thf", "solvent", E[5]); charge(pthf2, "volume", "1 mL", {"number": 1, "unit": "mL"}, E[5])
participant(r, "product12", "m12", "product", E[2] + E[3] + E[5], product_role="major product")
participant(r, "product13", "m13", "product", E[2] + E[3] + E[5], product_role="isolated side product")
pnh = participant(r, "nh4cl", "m_nh4cl", "workup_agent", E[5], stage="workup"); charge(pnh, "volume", "1 mL", {"number": 1, "unit": "mL"}, E[5])
pex = participant(r, "ether_extract", "m_ether", "workup_agent", E[5], stage="workup"); charge(pex, "volume", "3 × 20 mL", {"number": 20, "unit": "mL", "count": 3}, E[5])
pbr = participant(r, "brine", "m_brine", "workup_agent", E[5], stage="workup"); charge(pbr, "volume", "15 mL", {"number": 15, "unit": "mL"}, E[5])
participant(r, "na2so4", "m_na2so4", "workup_agent", E[5], stage="workup"); participant(r, "petrol_purification", "m_petrol", "purification_agent", E[5], stage="purification"); participant(r, "ether_purification", "m_ether", "purification_agent", E[5], stage="purification")
scale(r, "m4", "0.123 mmol methyleneaziridine 4", {"number": 0.123, "unit": "mmol"}, E[5]); solvents(r, "THF: 1 mL CuI suspension and 1 mL substrate solution; vinylmagnesium bromide supplied in THF", ["m_thf"], E[5])
step(r, "organocopper_preparation", "reaction", ["vinylmgbr", "cui", "thf_cui"], E[5], "add vinylmagnesium bromide dropwise to CuI suspension and wait", temperature=("−50 °C", {"number": -50, "unit": "°C"}), time=("10 min", {"number": 10, "unit": "min"}))
temp_conflict = fact(["Scheme 4: −50 °C → 0 °C, 1 h", "prose: warm to −20 °C over 1 h, then quench"], None, E[3] + E[5], status="conflict", alternatives=[fact("−50 °C → 0 °C", {"minimum": -50, "maximum": 0, "unit": "°C", "raw_text": "−50 °C → 0 °C"}, E[3]), fact("−50 °C to −20 °C", {"minimum": -50, "maximum": -20, "unit": "°C", "raw_text": "−50 °C to −20 °C"}, E[5])])
step(r, "substitution", "reaction", ["substrate4", "vinylmgbr", "cui", "thf_substrate"], E[3] + E[5], "add solution of 4 via cannula and warm over 1 h", temperature_fact=temp_conflict, time=("1 h", {"number": 1, "unit": "h"}))
step(r, "quench", "workup", ["nh4cl"], E[5], "quench with saturated aqueous NH4Cl, stir for 30 min, then warm to RT", temperature=("warm to RT", "r.t."), time=("30 min", {"number": 30, "unit": "min"}), atmosphere=False)
step(r, "extract_wash_dry", "workup", ["ether_extract", "brine", "na2so4"], E[5], "extract with ether, wash with brine, dry with Na2SO4, and concentrate", atmosphere=False)
step(r, "flash", "purification", ["petrol_purification", "ether_purification"], E[5], "flash chromatography (petrol/ether, 7:1 to 1:1)", atmosphere=False)
procedure(r, "Vinylmagnesium bromide was added dropwise to CuI in THF at −50 °C; after 10 min a solution of 4 in THF was added via cannula and the mixture was warmed to −20 °C over 1 h before quench.", E[5], addition="vinylmagnesium bromide to CuI; after 10 min add 4 solution via cannula")
outcome(r, "yield12", "yield", ["m12"], E[3] + E[5], "31%", {"number": 31, "unit": "%"}, "major product", "isolated yield", "isolated after flash chromatography")
outcome(r, "amount12", "product_amount", ["m12"], E[5], "7.0 mg", {"number": 7.0, "unit": "mg"}, "major product", basis="isolated material")
outcome(r, "yield13", "yield", ["m13"], E[3] + E[5], "23%", {"number": 23, "unit": "%"}, "isolated side product", "isolated yield", "isolated after flash chromatography")
outcome(r, "amount13", "product_amount", ["m13"], E[5], "6.0 mg", {"number": 6.0, "unit": "mg"}, "isolated side product", basis="isolated material")
outcome(r, "stereo12", "selectivity", ["m12"], E[2] + E[3], "new C-C bond cis to the deuterium atom", "cis-to-D; stereoretentive", "major product", basis="NMR experiments including diagnostic NOE", scope="relative disposition of newly formed C-C bond and deuterium")
note(r, "The azirine 13 may form by competing addition at the carbonyl followed by 1,4-vinylation of the resulting α,β-unsaturated ester.", E[2])
set_classification(r, E[2] + E[3] + E[5], [("C", "N", 1, 0, "broken", "exo vinylic C-N cleavage in the 12 branch"), ("C", "C", 0, 1, "formed", "new vinylic carbon-vinyl bond in 12"), ("N", "H", 0, 1, "hydrogen_gain", "oxazolidinone N-H in 12 after workup"), ("C", "C", 0, 1, "formed", "at least one new carbon-vinyl bond in the proposed 13 branch")],
                   "One experiment affords two products; the 12 branch is source-drawn, while complete mapping of the proposed 13 branch is unresolved.")
bundle["reactions"].append(r)

# Scheme 5 formal insertion/ligand coupling proposal.
r = reaction("scheme5_cu_insertion_coupling", "proposed_mechanism", E[3], "page 3, Scheme 5",
             [], "retentive ring-opening via formal C-N insertion and ligand coupling", "proposed organocopper cross-coupling mechanism",
             ["Intermediates 14 and 15 are proposed generic organocopper species, not isolated compounds.", "Counterion M and substituent R are generic and full oxidation-state/atom mapping is not established."])
participant(r, "substrate4", "m4", "reactant", E[3]); participant(r, "cuprate", "m_generic_cuprate", "reagent", E[3]); participant(r, "intermediate14", "m14", "product", E[3], product_role="proposed intermediate"); participant(r, "intermediate15", "m15", "product", E[3], product_role="proposed ligand-coupling product")
step(r, "formal_insertion", "reaction", ["substrate4", "cuprate", "intermediate14"], E[3], "formal Cu(I) insertion into the activated sp2 C-N bond", atmosphere=False)
step(r, "ligand_coupling", "reaction", ["intermediate14", "intermediate15"], E[3], "reductive elimination/ligand coupling", atmosphere=False)
procedure(r, "Retentive ring-opening via formal C-N insertion then ligand coupling from Cu(III) intermediate 14.", E[3])
set_classification(r, E[3], [("C", "N", 1, 0, "broken", "proposed activated C-N bond cleavage"), ("C", "C", 0, 1, "formed", "proposed ligand-coupling C-C bond")],
                   "Proposed generic organocopper pathway; no experimental atom mapping or isolated intermediates.")
bundle["reactions"].append(r)

# Non-synthetic source context: characterization, stereochemical analysis, and discussion are not extra experiments.
bundle["source_context"] = [
    {"context_id": f"{DOC}:general_information", "document_id": DOC, "kind": "general_experimental_information", "evidence_ids": E[3],
     "raw_text": "All reactions were carried out in oven-dried glassware and under an atmosphere of argon unless otherwise specified; analytical and chromatography conventions are reported."},
    {"context_id": f"{DOC}:fig1_noe_4", "document_id": DOC, "kind": "stereochemical_characterization", "compound_labels": ["4"], "evidence_ids": E[2],
     "raw_text": "Fig. 1 gives NMR chemical shifts and diagnostic NOE correlations for methyleneaziridine 4; 3JHH = 9.5 Hz supports the assigned relative conformation."},
    {"context_id": f"{DOC}:fig2_noe_11_12", "document_id": DOC, "kind": "stereochemical_characterization", "compound_labels": ["11", "12"], "evidence_ids": E[2] + E[3],
     "raw_text": "Fig. 2 gives diagnostic NOE correlations for 11 and 12, including important weak/absent correlations."},
    {"context_id": f"{DOC}:characterization_7_13", "document_id": DOC, "kind": "characterization_only", "compound_labels": ["7", "8", "8-hydroxycarbamate", "9", "4", "11", "12", "13"], "evidence_ids": E[4] + E[5],
     "raw_text": "Rf, IR, NMR, HRMS and melting-point data accompanying the preparative paragraphs are retained as characterization context and were not converted into separate synthesis records."},
    {"context_id": f"{DOC}:mechanistic_discussion", "document_id": DOC, "kind": "author_mechanistic_discussion", "evidence_ids": E[2] + E[3],
     "raw_text": "The authors conclude that the substitution is stereoretentive and rules out S_NVσ, while S_NVπ/carbocupration-fragmentation and formal Cu insertion/ligand coupling remain mechanistic possibilities."},
    {"context_id": f"{DOC}:reference_only_pages", "document_id": DOC, "kind": "references", "evidence_ids": E[5] + E[6],
     "raw_text": "The right column of page 5 and page 6 are references/footnotes; the small structure 16 on page 6 is a cited comparison compound, not a synthesis reported in this paper."},
]

visual_inventory = {
    "document_id": DOC,
    "pages_reviewed": [1, 2, 3, 4, 5, 6],
    "observations": [
        {"page": 1, "regions": ["title/abstract", "Scheme 1", "introduction/results prose"], "finding": "Scheme 1 is explicitly prior work; two generic literature arrows were separated from current experiments."},
        {"page": 2, "regions": ["Scheme 2", "Fig. 1", "Scheme 3", "results prose"], "finding": "Two hypothetical S_NV modes, stereochemical NMR context, and five Scheme 3 transformations were inventoried."},
        {"page": 3, "regions": ["Scheme 4", "Fig. 2", "Scheme 5", "general information", "start of compound 7 procedure"], "finding": "Two current S_NV experiments, one proposed Cu mechanism, and the beginning of the deuteration procedure were identified."},
        {"page": 4, "regions": ["compound 7 continuation", "procedures for 8, hydroxycarbamate/9, 4, and 11"], "finding": "Preparative operations, quantities, yields, workups, and characterization boundaries were reviewed across both columns."},
        {"page": 5, "regions": ["11 characterization continuation", "12/13 procedure", "references"], "finding": "The shared 12/13 experiment and two isolated yields were captured; peak lists remain characterization context."},
        {"page": 6, "regions": ["references and footnotes 18–23", "comparison structure 16"], "finding": "No current-paper synthetic experiment; structure 16 is cited comparison data only."},
    ],
    "semantic_accuracy_verified": False,
}

stereo_diagnostics = {}
for molecule in bundle["molecules"]:
    smi = molecule["smiles"]["normalized_value"]
    if smi:
        stereo_diagnostics[molecule["molecule_id"]] = st.stereo_inspect(smi)

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/visual_inventory_pages_0001_0006.json").write_text(json.dumps(visual_inventory, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/stereo_diagnostics.json").write_text(json.dumps(stereo_diagnostics, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "source_context": len(bundle["source_context"])}, ensure_ascii=False))
