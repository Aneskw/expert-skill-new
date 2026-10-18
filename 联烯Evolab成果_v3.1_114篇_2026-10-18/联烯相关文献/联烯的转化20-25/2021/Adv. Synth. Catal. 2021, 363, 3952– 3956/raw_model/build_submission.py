#!/usr/bin/env python3
"""Build the source-specific v3.1 bundle for doc_e6c9ddf245b60941.

The code is intentionally source-specific and reproducible.  It does not read old
reaction records; it starts from the standard-init template and the supplied page
evidence.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_e6c9ddf245b60941"
bundle = json.loads((RUN / "submitted.json").read_text())
assert bundle["schema_version"] == "3.1"
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []


def ev(page, image=True):
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"] if image else [f"{DOC}:page{page}:text"]


def F(raw=None, norm=None, evidence=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, norm, evidence, method=method, derivation=derivation,
                   status=status)


def D(norm, evidence, derivation):
    return F(None, norm, evidence, "structure_derivation", derivation)


def Q(number, unit, raw=None):
    return {"number": number, "unit": unit, **({"raw_text": raw} if raw else {})}


def set_statuses(m, chirality, geometry, representation, evidence, note):
    s = m["stereochemistry"]
    s["chirality_status"] = D(chirality, evidence, note)
    s["geometry_status"] = D(geometry, evidence, note)
    s["representation_status"] = D(representation, evidence, note)


def molecule(mid, label, name, evidence, smiles=None, chirality="unspecified",
             geometry="not_applicable", representation="partial", stereo_note=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(name, name, evidence)
    m["names"]["english"] = F(name, name, evidence)
    m["cas_numbers"] = F()
    m["identity_note"] = F(
        "Candidate identity transcribed from the supplied paper; no external identity lookup was performed.",
        "Candidate identity transcribed from source; no external lookup.", evidence,
    )
    m["external_properties"] = []
    note = stereo_note or "Source drawing contains no explicit stereochemical descriptor; status is conservative."
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = D(smiles, evidence, "Candidate SMILES transcribed from the source drawing/name; RDKit parsing is not source verification.")
        m["stereochemistry"]["isomeric_smiles"] = D(diag["canonical_isomeric_smiles"], evidence, "RDKit canonicalization of the candidate source-derived graph.")
        m["stereochemistry"]["identity_key"] = D(diag["stereo_graph_key"], evidence, "Stereo graph key calculated from the candidate SMILES; not a sample identity.")
    set_statuses(m, chirality, geometry, representation, evidence, note)
    return m


def generic_molecule(mid, label, name, evidence, fixed, variables):
    m = molecule(mid, label, name, evidence, smiles=None, chirality="unresolved",
                 geometry="unresolved", representation="not_encoded",
                 stereo_note="Generic R-group drawing does not define a single stereochemical identity.")
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": F(fixed, fixed, evidence),
        "variables": [
            {"label": v, "attachment_points": F(f"source-drawn attachment for {v}", f"source-drawn attachment for {v}", evidence),
             "definition": F(f"Variable {v} as printed", f"Variable {v} as printed", evidence),
             "allowed_values": F(), "constraints": F()}
            for v in variables
        ],
        "instance_mapping": [],
        "reason": F("R-group/general structure in source", "R-group/general structure in source", evidence),
    }
    return m


M = {}
def add(m):
    M[m["molecule_id"]] = m


# Current-work substrates and products.  Terminal-CH2 allenes are not assigned an
# allene configuration; their carbinol center remains unspecified.
allenes = {
    "1a": "C=C=C(C)C(O)c1ccccc1",
    "1b": "C=C=C(C)C(O)c1ccc(OC)cc1",
    "1c": "C=C=C(C)C(O)c1ccc(C)cc1",
    "1d": "C=C=C(C)C(O)c1ccc(C(=O)OC)cc1",
    "1e": "C=C=C(C)C(O)c1ccccc1C",
    "1f": "C=C=C(C)C(O)c1cccc(OC)c1",
    "1g": "C=C=C(C)C(O)c1ccc(OC)cc1OC",
    "1h": "C=C=C(C)C(O)c1cccs1",
    "1i": "C=C=C(C)C(O)c1ccsc1",
    "1j": "C=C=C(c1ccccc1)C(O)c1ccc(OC)cc1",
    "1k": "C=C=C(c1ccccc1)C(O)c1ccc(C)cc1",
    "1l": "C=C=C(c1ccccc1)C(O)c1ccc(Br)cc1",
    "1m": "CCCCCC(O)C=C=Cc1ccccc1",
    "1n": "CC(O)(c1ccccc1)C=C=C(C)C",
}
for label, smi in allenes.items():
    page = 2 if label == "1a" else 3
    name = {"1m": "1-phenylnona-1,2-dien-4-ol", "1n": "5-methyl-2-phenylhexa-3,4-dien-2-ol"}.get(label, f"allenol {label}")
    axial = "The 1,3-disubstituted allene and carbinol center are drawn without configuration." if label == "1m" else "The source gives no configuration; terminal-CH2 or symmetrically substituted allene axis is not assigned, and any carbinol center remains unspecified."
    add(molecule(f"{DOC}:mol_{label}", label, name, ev(page), smi,
                 chirality="unspecified", geometry="not_applicable",
                 representation="partial", stereo_note=axial))

products = {
    "3a": "CC1=C(S(=O)(=O)c2ccccc2)COC1c1ccccc1",
    "3b": "CC1=C(S(=O)(=O)c2ccccc2)COC1c1ccc(OC)cc1",
    "3c": "CC1=C(S(=O)(=O)c2ccccc2)COC1c1ccc(C)cc1",
    "3d": "CC1=C(S(=O)(=O)c2ccccc2)COC1c1ccc(C(=O)OC)cc1",
    "3e": "CC1=C(S(=O)(=O)c2ccccc2)COC1c1ccccc1C",
    "3f": "CC1=C(S(=O)(=O)c2ccccc2)COC1c1cccc(OC)c1",
    "3g": "CC1=C(S(=O)(=O)c2ccccc2)COC1c1ccc(OC)cc1OC",
    "3h": "CC1=C(S(=O)(=O)c2ccccc2)COC1c1cccs1",
    "3i": "CC1=C(S(=O)(=O)c2ccccc2)COC1c1ccsc1",
    "3j": "c1ccc(cc1)C1=C(S(=O)(=O)c2ccccc2)COC1c1ccc(OC)cc1",
    "3k": "c1ccc(cc1)C1=C(S(=O)(=O)c2ccccc2)COC1c1ccc(C)cc1",
    "3l": "c1ccc(cc1)C1=C(S(=O)(=O)c2ccccc2)COC1c1ccc(Br)cc1",
    "3m": "CC1=C(S(=O)(=O)c2ccc(C)cc2)COC1c1ccc(OC)cc1",
    "3n": "CC1=C(S(=O)(=O)c2ccc(C)cc2)COC1c1ccc(Br)cc1",
    "3o": "c1ccc(cc1)C1=C(S(=O)(=O)c2ccc(C)cc2)COC1c1ccc(OC)cc1OC",
    "3p": "CC1=C(S(=O)(=O)c2ccc(C)cc2)COC1c1cccs1",
}
for label, smi in products.items():
    page = 2 if label == "3a" else 3
    add(molecule(f"{DOC}:mol_{label}", label, f"4-(arylsulfonyl)-2,5-dihydrofuran {label}", ev(page), smi,
                 chirality="unspecified", geometry="not_applicable", representation="partial",
                 stereo_note="The ring carbon bearing aryl/H is potentially stereogenic; no wedge, R/S, ee, er, dr, or racemic statement is printed."))
add(molecule(f"{DOC}:mol_4a", "4a", "fully aromatic furan 4a", ev(3),
             "c1ccc(cc1)c1c(S(=O)(=O)c2ccccc2)coc1c1ccc(OC)cc1",
             chirality="achiral", geometry="not_applicable", representation="encoded",
             stereo_note="Candidate aromatic graph has no assigned stereogenic center or axis."))

# Sulfinates, catalysts, oxidants, solvents, and procedure agents.
simple = {
    "2a": ("2a", "sodium benzenesulfinate", "[Na+].[O-]S(=O)c1ccccc1"),
    "2b": ("2b", "sodium p-toluenesulfinate", "[Na+].[O-]S(=O)c1ccc(C)cc1"),
    "CuOAc": ("CuOAc", "copper(I) acetate", "[Cu+].CC(=O)[O-]"),
    "CuOAc2": ("Cu(OAc)2", "copper(II) acetate", "[Cu+2].CC(=O)[O-].CC(=O)[O-]"),
    "CuCl2": ("CuCl2", "copper(II) chloride", "[Cu+2].[Cl-].[Cl-]"),
    "FeCl2": ("FeCl2", "iron(II) chloride", "[Fe+2].[Cl-].[Cl-]"),
    "AgNO3": ("AgNO3", "silver nitrate", "[Ag+].[O-][N+](=O)[O-]"),
    "Ag2CO3": ("Ag2CO3", "silver carbonate", "[Ag+].[Ag+].[O-]C([O-])=O"),
    "AgOAc": ("AgOAc", "silver acetate", "[Ag+].CC(=O)[O-]"),
    "K2S2O8": ("K2S2O8", "potassium persulfate", "[K+].[K+].[O-]S(=O)(=O)OOS(=O)(=O)[O-]"),
    "TBPB": ("TBPB", "tert-butyl peroxybenzoate", "CC(C)(C)OOC(=O)c1ccccc1"),
    "O2": ("O2", "molecular oxygen", "O=O"),
    "TEMPO": ("TEMPO", "2,2,6,6-tetramethylpiperidin-1-oxyl", "CC1(C)CCCC(C)(C)N1[O]"),
    "MeCN": ("MeCN", "acetonitrile", "CC#N"),
    "DMF": ("DMF", "N,N-dimethylformamide", "CN(C)C=O"),
    "DCE": ("DCE", "1,2-dichloroethane", "ClCCCl"),
    "THF": ("THF", "tetrahydrofuran", "C1CCOC1"),
    "EtOAc": ("ethyl acetate", "ethyl acetate", "CCOC(C)=O"),
    "MgSO4": ("MgSO4", "magnesium sulfate", "[Mg+2].[O-]S(=O)(=O)[O-]"),
    "hexane": ("n-hexane", "n-hexane", "CCCCCC"),
}
for key, (label, name, smi) in simple.items():
    evidence = ev(3) if key in ("2b", "TEMPO") else ev(4) if key in ("EtOAc", "MgSO4", "hexane") else ev(2)
    simple_chirality = "unspecified" if key in ("2a", "2b") else "achiral"
    simple_representation = "partial" if key in ("2a", "2b") else "encoded"
    add(molecule(f"{DOC}:mol_{key}", label, name, evidence, smi,
                 chirality=simple_chirality, geometry="not_applicable", representation=simple_representation,
                 stereo_note="For sulfinate sulfur, the source does not assign configuration; other simple candidate reagent/solvent graphs have no assigned stereogenic element."))
add(generic_molecule(f"{DOC}:mol_brine", "brine", "brine", ev(4), "aqueous sodium chloride mixture", []))

# General and literature-only structures from Scheme 1.
add(generic_molecule(f"{DOC}:mol_generic_allenol", "1", "generic α-allenol 1", ev(2), "R1/R2-substituted α-allenol", ["R1", "R2"]))
add(generic_molecule(f"{DOC}:mol_generic_sulfinate", "2", "generic sodium arenesulfinate 2", ev(2), "ArSO2Na", ["Ar"]))
add(generic_molecule(f"{DOC}:mol_generic_product", "3", "generic 4-(arylsulfonyl)-2,5-dihydrofuran 3", ev(2), "R1/R2/Ar-substituted dihydrofuran sulfone", ["R1", "R2", "Ar"]))
for n, fixed, variables in [
    ("lit_a_sub", "allene-alkyne alcohol bearing R1/R2", ["R1", "R2"]),
    ("lit_a_prod", "benzannulated sulfonylated product", ["R1", "R2", "Ar"]),
    ("lit_b_sub", "terminal allenol", []),
    ("lit_b_prod", "hydrosulfinated/arylated alcohol", ["Ar", "ArI"]),
    ("lit_c_sub", "R1/R2-substituted α-allenol", ["R1", "R2"]),
    ("lit_c_prod", "sulfonylation-semipinacol aldehyde", ["R1", "R2", "Ar"]),
    ("lit_e_sub", "sulfonylallenol", ["R1", "R2", "R3"]),
    ("lit_e_prod", "sulfonyl-dihydrofuran", ["R1", "R2", "R3"]),
]:
    add(generic_molecule(f"{DOC}:mol_{n}", n, n.replace("_", " "), ev(2), fixed, variables))


def participant(pid, mid, role, evidence, stage="reaction", mass=None, amount=None,
                volume=None, equiv=None, molpct=None, ref=None, product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage)
    if product_role is not None:
        p["product_role"] = F(product_role, product_role, evidence)
    for key, val in (("mass", mass), ("amount", amount), ("volume", volume),
                     ("equivalents", equiv), ("mol_percent", molpct)):
        if val is not None:
            p["charges"][key] = F(val[0], val[1], evidence)
    if ref:
        p["equivalence_reference"] = F(ref, ref, evidence)
    return p


def base_reaction(suffix, evidence, record_type, categories, location, name, transformation):
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories, categories, evidence)
    r["location"] = F(location, location, evidence)
    r["reaction_name"] = F(name, name, evidence)
    r["transformation_type"] = F(None, transformation, evidence, "inference",
                                  "Normalized transformation description based on the source caption/prose.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def set_classification(r, kind, evidence):
    c = r["classification"]
    c["mapped_reaction_smiles"] = F()
    c["model_votes"] = []
    if kind == "cascade":
        changes = [
            {"change_id": "b1", "atoms": [{"element": "C"}, {"element": "O"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"element": "C"}, {"element": "S"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b3", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
            {"change_id": "b4", "atoms": [{"element": "O"}, {"element": "H"}], "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "hydrogen_count": 1},
        ]
        labels = [{"category_number": 3, "code": "C-O", "supporting_change_ids": ["b1"]},
                  {"category_number": 4, "code": "C-S", "supporting_change_ids": ["b2"]}]
        deriv = "Partial class-level changes inferred from the source conclusion and drawings; exact atom mapping and coproducts are unresolved."
        c["bond_changes"] = [D(x, evidence, deriv) for x in changes]
        c["labels"] = D(labels, evidence, deriv)
        c["coverage"] = D("partial", evidence, deriv)
        c["limitations"] = ["Exact atom identities and complete atom mapping were not asserted.", "Hydrogen/coproduct accounting is incomplete."]
    elif kind == "aromatization":
        changes = [
            {"change_id": "b1", "atoms": [{"element": "C"}, {"element": "O"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"element": "C"}, {"element": "S"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b3", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 1, "after_order": 2, "change_type": "order_increased"},
            {"change_id": "b4", "atoms": [{"element": "C"}, {"element": "H"}], "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "hydrogen_count": 1},
        ]
        labels = [{"category_number": 3, "code": "C-O", "supporting_change_ids": ["b1"]},
                  {"category_number": 4, "code": "C-S", "supporting_change_ids": ["b2"]},
                  {"category_number": 13, "code": "DB-C-C", "supporting_change_ids": ["b3"]}]
        deriv = "Class-level cyclization/sulfonylation plus aromatization inferred from Scheme 3b; exact mapping is unresolved."
        c["bond_changes"] = [D(x, evidence, deriv) for x in changes]
        c["labels"] = D(labels, evidence, deriv)
        c["coverage"] = D("partial", evidence, deriv)
        c["limitations"] = ["Exact atom mapping, oxidation partner, and hydrogen balance are unresolved."]
    elif kind == "not_applicable":
        deriv = "No characterized target transformation was observed; intended bond formation is not assigned."
        c["bond_changes"] = []
        c["labels"] = D([], evidence, deriv)
        c["coverage"] = D("not_applicable", evidence, deriv)
        c["limitations"] = ["Negative or inhibited experiment; no intended transformation labels inferred."]
    else:
        deriv = "Source drawing is generic, cited prior work, or gives an unidentified mixture; no complete classification asserted."
        c["bond_changes"] = []
        c["labels"] = D([], evidence, deriv)
        c["coverage"] = D("unresolved", evidence, deriv)
        c["limitations"] = ["No safe exact atom mapping/classification from the supplied generic or unidentified structures."]


def reaction_step(r, evidence, pids, temp=None, time=None, operation="reaction mixture held under the stated conditions", stage="reaction", atmosphere=None):
    s = st.new_step(f"step_{len(r['steps'])+1}", evidence, stage)
    s["participant_ids"] = list(pids)
    s["operation"] = F(operation, operation, evidence)
    if temp is not None: s["temperature"] = F(temp[0], temp[1], evidence)
    if time is not None: s["time"] = F(time[0], time[1], evidence)
    if atmosphere is not None: s["atmosphere"] = F(atmosphere, atmosphere, evidence)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, evidence, raw, norm, basis=None, yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = F(raw, norm, evidence)
    if basis is not None: o["measurement_basis"] = F(basis, basis, evidence)
    if yield_type is not None: o["yield_type"] = F(yield_type, yield_type, evidence)
    if scope is not None: o["scope"] = F(scope, scope, evidence)
    if kind in ("yield", "product_amount"):
        o["product_role"] = F("product", "product", evidence)
    r["outcomes"].append(o)


R = []
def add_r(r): R.append(r)

# Scheme 1: four cited examples and one current general scheme.
lit_specs = [
    ("scheme1a_literature_jiang2017", "lit_a_sub", "lit_a_prod", "Benzannulation-sulfonylation", "Jiang, 2017; IBX then ArN2BF4/[DABSO]"),
    ("scheme1b_literature_gong2018", "lit_b_sub", "lit_b_prod", "Hydrosulfination-arylation", "Gong, 2018; TsNHNH2 + ArI, [Pd]"),
    ("scheme1c_literature_almendros2020", "lit_c_sub", "lit_c_prod", "Sulfonylation-semipinacol", "Almendros, 2020; ArN2BF4 + DABSO, Ru(II), LED lamp"),
    ("scheme1e_literature_harmata2016", "lit_e_sub", "lit_e_prod", "Cycloisomerization of sulfonylallenols", "Harmata, 2016; AgF"),
]
for suffix, sub, prod, name, note in lit_specs:
    e = ev(2)
    r = base_reaction(suffix, e, "literature_example", ["other"], "Scheme 1 cited prior-work arrow", name, name.lower())
    r["participants"] = [participant("substrate", f"{DOC}:mol_{sub}", "reactant", e), participant("product", f"{DOC}:mol_{prod}", "product", e, product_role="reported prior-work product")]
    reaction_step(r, e, ["substrate", "product"], operation=note)
    r["notes"].append(F(note, note, e))
    r["review"]["unresolved"] = ["Scheme 1 is a cited literature summary; experimental quantities and complete structures are not provided on this page."]
    set_classification(r, "unresolved", e); add_r(r)

e = ev(2)
general = base_reaction("scheme1d_current_generic", e, "generic_scheme", ["other"], "Scheme 1d; repeated headers in Schemes 2-3", "Cycloetherification-sulfonylation", "cascade cycloetherification/sulfonylation")
general["participants"] = [participant("allenol", f"{DOC}:mol_generic_allenol", "reactant", e), participant("sulfinate", f"{DOC}:mol_generic_sulfinate", "reactant", e), participant("product", f"{DOC}:mol_generic_product", "product", e, product_role="desired product"), participant("cu", f"{DOC}:mol_CuOAc2", "catalyst", e), participant("ag", f"{DOC}:mol_AgNO3", "reagent", e)]
general["review"]["unresolved"] = ["R-group scheme is not one concrete experiment; variable structures and quantities are not fully specified here."]
set_classification(general, "cascade", e); add_r(general)
GEN_ID = general["reaction_id"]


def table_reaction(entry, catalyst, oxidant, solvent, temp, hours, result):
    e = ev(2)
    suffix = f"table1_entry{entry}"
    cats = ["optimization"] + (["substrate_scope"] if entry == 2 else [])
    r = base_reaction(suffix, e, "specific_experiment", cats, f"Table 1 entry {entry}", "Cycloetherification/sulfonylation of allenol 1a", "cascade cycloetherification/sulfonylation")
    r["general_scheme_ids"] = [GEN_ID]
    mids = {"s": f"{DOC}:mol_1a", "sf": f"{DOC}:mol_2a", "p": f"{DOC}:mol_3a", "cat": f"{DOC}:mol_{catalyst}", "ox": f"{DOC}:mol_{oxidant}", "solv": f"{DOC}:mol_{solvent}"}
    r["participants"] = [
        participant("allenol_1a", mids["s"], "reactant", e, amount=("0.15 mmol", Q(0.15, "mmol")), equiv=("1 equiv", Q(1, "equiv")), ref=mids["s"]),
        participant("sulfinate_2a", mids["sf"], "reactant", e, amount=("0.30 mmol", Q(0.30, "mmol")), equiv=("2 equiv", Q(2, "equiv")), ref=mids["s"]),
        participant("product_3a", mids["p"], "product", e, product_role="desired product"),
        participant("catalyst", mids["cat"], "catalyst", e, molpct=("10 mol%", Q(10, "mol%")), ref=mids["s"]),
        participant("oxidant", mids["ox"], "reagent", e, molpct=("200 mol%", Q(200, "mol%")), ref=mids["s"]),
        participant("solvent", mids["solv"], "solvent", e, volume=("1 mL", Q(1, "mL"))),
    ]
    r["reaction_scale"] = F("1a (0.15 mmol)", {"reference_molecule_id": mids["s"], "amount": Q(0.15, "mmol")}, e)
    r["solvent_conditions"]["composition"] = F(solvent, [solvent], e)
    r["solvent_conditions"]["total_volume"] = F("1 mL", Q(1, "mL"), e)
    reaction_step(r, e, [p["participant_id"] for p in r["participants"]], (f"{temp} °C", Q(temp, "°C")), (f"{hours} h", Q(hours, "h")), atmosphere=None)
    r["special_conditions"]["reactor"] = F("sealed tube", "sealed tube", e)
    r["procedure"]["original_text"] = F("Table 1 footnote: sealed tube using 1a (0.15 mmol), 2a (0.30 mmol), oxidant (200 mol%) and catalyst (10 mol%) in solvent (1 mL); progress followed by TLC.", "Table 1 footnote procedure", e)
    if result == "trace":
        outcome(r, "trace_3a", "trace", [mids["p"]], e, "trace", "trace", basis="pure, isolated product with correct analytical and spectral data", scope={"table_column": "yield (%)"})
    else:
        outcome(r, "yield_3a", "yield", [mids["p"]], e, f"{result}%", Q(result, "%"), basis="pure, isolated product with correct analytical and spectral data", yield_type="isolated", scope={"basis_substrate": "1a"})
    r["review"]["unresolved"] = ["Atmosphere and pressure are not reported; complete workup is absent from Table 1."]
    set_classification(r, "cascade", e)
    return r


table_rows = [
    (1,"CuOAc","AgNO3","MeCN",100,3,45),(2,"CuOAc2","AgNO3","MeCN",100,2,68),
    (3,"CuCl2","AgNO3","MeCN",100,3,37),(4,"FeCl2","AgNO3","MeCN",100,3,"trace"),
    (5,"CuOAc2","Ag2CO3","MeCN",100,2,49),(6,"CuOAc2","AgOAc","MeCN",100,2,46),
    (7,"CuOAc2","K2S2O8","MeCN",100,3,"trace"),(8,"CuOAc2","AgNO3","DMF",100,2.5,51),
    (9,"CuOAc2","AgNO3","DCE",100,3,"trace"),(10,"CuOAc2","AgNO3","THF",100,3,21),
    (11,"CuOAc2","AgNO3","MeCN",40,3,"trace"),(12,"CuOAc2","TBPB","MeCN",100,3,18),
    (13,"CuOAc2","O2","MeCN",100,3,"trace"),
]
for row in table_rows: add_r(table_reaction(*row))


def scope_reaction(label, substrate, sulfinate, product, percent, hours, detailed=False):
    e = ev(3) + (ev(4) if detailed else [])
    r = base_reaction(f"scheme_scope_{label}", e, "specific_experiment", ["substrate_scope"], f"Scheme {'2' if label <= '3i' else '3a'} product {label}", "Cooperative copper-silver catalyzed cycloetherification/sulfonylation", "cascade cycloetherification/sulfonylation")
    r["general_scheme_ids"] = [GEN_ID]
    smid, sfmid, pmid = f"{DOC}:mol_{substrate}", f"{DOC}:mol_{sulfinate}", f"{DOC}:mol_{product}"
    r["participants"] = [
        participant("substrate", smid, "reactant", e),
        participant("sulfinate", sfmid, "reactant", e, equiv=("2 equiv", Q(2, "equiv")), ref=smid),
        participant("product", pmid, "product", e, product_role="isolated product"),
        participant("catalyst", f"{DOC}:mol_CuOAc2", "catalyst", e, molpct=("10 mol%", Q(10, "mol%")), ref=smid),
        participant("oxidant", f"{DOC}:mol_AgNO3", "reagent", e, equiv=("2 equiv", Q(2, "equiv")), ref=smid),
        participant("solvent", f"{DOC}:mol_MeCN", "solvent", e),
    ]
    r["solvent_conditions"]["composition"] = F("MeCN", ["MeCN"], e)
    reaction_step(r, e, [p["participant_id"] for p in r["participants"]], ("100 °C", Q(100,"°C")), (f"{hours} h", Q(hours,"h")))
    r["special_conditions"]["reactor"] = F("sealed tube", "sealed tube", e)
    outcome(r, f"yield_{product}", "yield", [pmid], e, f"{percent}%", Q(percent,"%"), basis="isolated yield", yield_type="isolated")
    if detailed:
        # Replace generic charges with the page-4 concrete procedure for 3b.
        vals = {"substrate": (("34.8 mg",Q(34.8,"mg")),("0.18 mmol",Q(0.18,"mmol"))),
                "sulfinate": (("59.1 mg",Q(59.1,"mg")),("0.36 mmol",Q(0.36,"mmol"))),
                "catalyst": (("3.3 mg",Q(3.3,"mg")),("0.018 mmol",Q(0.018,"mmol"))),
                "oxidant": (("61.2 mg",Q(61.2,"mg")),("0.36 mmol",Q(0.36,"mmol")))}
        for p in r["participants"]:
            if p["participant_id"] in vals:
                p["charges"]["mass"] = F(vals[p["participant_id"]][0][0], vals[p["participant_id"]][0][1], ev(4))
                p["charges"]["amount"] = F(vals[p["participant_id"]][1][0], vals[p["participant_id"]][1][1], ev(4))
            if p["participant_id"] == "solvent": p["charges"]["volume"] = F("2 mL",Q(2,"mL"),ev(4))
        r["reaction_scale"] = F("allenol 1b (0.18 mmol)", {"reference_molecule_id": smid, "amount":Q(0.18,"mmol")}, ev(4))
        r["solvent_conditions"]["total_volume"] = F("2 mL",Q(2,"mL"),ev(4))
        r["procedure"]["original_text"] = F("Typical Procedure for synthesis of 3b, from reagent addition through chromatography and characterization.", "page-4 typical procedure for 3b", ev(4))
        reaction_step(r, ev(4), ["product"], operation="cool to room temperature", stage="workup")
        for pid, mid, role, raw, norm in [
            ("ethyl_acetate",f"{DOC}:mol_EtOAc","workup_agent","ethyl acetate (3 × 2 mL)",{"count":3,"volume_each":Q(2,"mL")}),
            ("brine",f"{DOC}:mol_brine","workup_agent","brine","brine"),
            ("drying_agent",f"{DOC}:mol_MgSO4","workup_agent","anhydrous MgSO4","anhydrous MgSO4"),
            ("hexane",f"{DOC}:mol_hexane","purification_agent","n-hexane","n-hexane"),
        ]:
            p = participant(pid,mid,role,ev(4),stage="workup" if role=="workup_agent" else "purification")
            if pid=="ethyl_acetate": p["charges"]["volume"] = F(raw,Q(6,"mL"),ev(4),"inference","Total extraction volume calculated as 3 × 2 mL; the three portions remain in raw_value.")
            r["participants"].append(p)
        # Ethyl acetate is also the chromatography co-eluent; one participant is reused in purification.
        reaction_step(r, ev(4), ["ethyl_acetate"], operation="extract with ethyl acetate (3 × 2 mL); separate organic layer", stage="workup")
        reaction_step(r, ev(4), ["brine"], operation="wash organic extract with brine", stage="workup")
        reaction_step(r, ev(4), ["drying_agent"], operation="dry over anhydrous MgSO4 and concentrate under reduced pressure", stage="workup")
        reaction_step(r, ev(4), ["hexane","ethyl_acetate"], operation="chromatography, n-hexane/ethyl acetate 7:1", stage="purification")
        outcome(r,"amount_3b","product_amount",[pmid],ev(4),"37 mg",Q(37,"mg"),basis="after chromatography")
        r["notes"].append(F("Scheme 2 reports 3 h; the Typical Procedure says until starting material disappeared by TLC without repeating a numerical duration.", "Scheme time and prose endpoint both retained", e))
    r["review"]["unresolved"] = ["Pressure and atmosphere are unreported."] + ([] if detailed else ["Substrate scale, solvent volume, workup, and product mass are not reported on the supplied pages."])
    set_classification(r,"cascade",e); return r


# 3a is already Table 1 entry 2; add its Scheme 2 evidence instead of duplicating.
entry2 = next(r for r in R if r["reaction_id"].endswith("table1_entry2"))
entry2["evidence_ids"] += [x for x in ev(3) if x not in entry2["evidence_ids"]]
entry2["location"] = F(["Table 1 entry 2","Scheme 2 product 3a"], ["Table 1 entry 2","Scheme 2 product 3a"], ev(2)+ev(3))
entry2["notes"].append(F("Table 1 entry 2 and Scheme 2 product 3a both report 68% at 100 °C for 2 h under the same catalyst/oxidant/solvent system; treated as one experiment.", "cross-page duplicate merged", ev(2)+ev(3)))
for label, pct, hrs in [("3b",62,3),("3c",59,2),("3d",63,1),("3e",66,2),("3f",63,2.5),("3g",69,16),("3h",68,16),("3i",61,1.5)]:
    add_r(scope_reaction(label,"1"+label[1:],"2a",label,pct,hrs,detailed=(label=="3b")))
for label, sub, sulf, pct, hrs in [
    ("3j","1j","2a",71,3),("3k","1k","2a",50,2),("3l","1l","2a",68,1.5),
    ("3m","1j","2b",59,1.5),("3n","1l","2b",59,1.5),("3o","1j","2b",65,2.5),("3p","1h","2b",74,2),
]: add_r(scope_reaction(label,sub,sulf,label,pct,hrs))

# Nonproductive 1m/1n experiments under the Scheme 3a standard conditions.
for sub in ("1m","1n"):
    e=ev(3); r=base_reaction(f"scheme3_nonproductive_{sub}",e,"specific_experiment",["substrate_scope"],f"Scheme 3a non-productive precursor {sub}","Cooperative copper-silver catalyzed cycloetherification/sulfonylation attempt","unproductive scope attempt")
    smid=f"{DOC}:mol_{sub}"; sfmid=f"{DOC}:mol_2a"
    r["general_scheme_ids"]=[GEN_ID]
    r["participants"]=[participant("substrate",smid,"reactant",e),participant("sulfinate",sfmid,"reactant",e,equiv=("2 equiv",Q(2,"equiv")),ref=smid),participant("catalyst",f"{DOC}:mol_CuOAc2","catalyst",e,molpct=("10 mol%",Q(10,"mol%")),ref=smid),participant("oxidant",f"{DOC}:mol_AgNO3","reagent",e,equiv=("2 equiv",Q(2,"equiv")),ref=smid),participant("solvent",f"{DOC}:mol_MeCN","solvent",e)]
    r["solvent_conditions"]["composition"]=F("MeCN",["MeCN"],e)
    reaction_step(r,e,[p["participant_id"] for p in r["participants"]],("100 °C",Q(100,"°C")),operation="exposed to the Scheme 3a standard conditions in a sealed tube")
    r["special_conditions"]["reactor"]=F("sealed tube","sealed tube",e)
    outcome(r,"nonproductive_result","other",[],e,"non-productive precursor; reaction was messy and did not produce the target sulfonyl-dihydrofuran","target sulfonyl-dihydrofuran not produced; messy reaction")
    r["review"]["unresolved"]=["Reaction time, conversion, recovered substrate, byproducts, and detection/measurement basis are unreported."]
    set_classification(r,"not_applicable",e); add_r(r)

# Prolonged-time aromatic furan 4a.
e=ev(3); r=base_reaction("scheme3b_furan_4a",e,"specific_experiment",["substrate_scope","other"],"Scheme 3b","Formation of fully aromatic furan 4a","cycloetherification/sulfonylation with aromatization")
smid=f"{DOC}:mol_1j"; pmid=f"{DOC}:mol_4a"
r["general_scheme_ids"]=[GEN_ID]
r["participants"]=[participant("substrate",smid,"reactant",e),participant("sulfinate",f"{DOC}:mol_2a","reactant",e,equiv=("2 equiv",Q(2,"equiv")),ref=smid),participant("product",pmid,"product",e,product_role="fully aromatic furan"),participant("catalyst",f"{DOC}:mol_CuOAc2","catalyst",e,molpct=("10 mol%",Q(10,"mol%")),ref=smid),participant("oxidant",f"{DOC}:mol_AgNO3","reagent",e,equiv=("2 equiv",Q(2,"equiv")),ref=smid),participant("solvent",f"{DOC}:mol_MeCN","solvent",e)]
r["solvent_conditions"]["composition"]=F("MeCN",["MeCN"],e)
reaction_step(r,e,[p["participant_id"] for p in r["participants"]],("100 °C",Q(100,"°C")),("15 h",Q(15,"h")))
r["special_conditions"]["reactor"]=F("sealed tube","sealed tube",e)
outcome(r,"yield_4a","yield",[pmid],e,"60%",Q(60,"%"),basis="isolated yield",yield_type="isolated")
r["review"]["unresolved"]=["Scale, solvent volume, atmosphere, product mass, and workup are unreported."]
set_classification(r,"aromatization",e); add_r(r)

# Scheme 4 controls.
controls=[
    ("control_without_silver","without silver salt",True,False,None,"complex mixture"),
    ("control_without_copper","without copper catalyst",False,True,None,"cycloisomerization was detected"),
    ("control_tempo","TEMPO (2 equiv)",True,True,"TEMPO","reaction completely inhibited"),
]
for suffix,desc,has_cu,has_ag,extra,result_note in controls:
    e=ev(3); r=base_reaction(suffix,e,"specific_experiment",["control","mechanistic_study"],f"Scheme 4: {desc}","Control experiment for radical/copper-silver pathway","mechanistic control")
    smid=f"{DOC}:mol_1a"; pmid=f"{DOC}:mol_3a"
    parts=[participant("allenol",smid,"reactant",e),participant("sulfinate",f"{DOC}:mol_2a","reactant",e,equiv=("2 equiv",Q(2,"equiv")),ref=smid),participant("product",pmid,"product",e,product_role="desired product"),participant("solvent",f"{DOC}:mol_MeCN","solvent",e)]
    if has_cu: parts.append(participant("catalyst",f"{DOC}:mol_CuOAc2","catalyst",e,molpct=("10 mol%",Q(10,"mol%")),ref=smid))
    if has_ag: parts.append(participant("oxidant",f"{DOC}:mol_AgNO3","reagent",e,equiv=("2 equiv",Q(2,"equiv")),ref=smid))
    if extra: parts.append(participant("tempo",f"{DOC}:mol_TEMPO","reagent",e,equiv=("2 equiv",Q(2,"equiv")),ref=smid))
    r["participants"]=parts;r["solvent_conditions"]["composition"]=F("MeCN",["MeCN"],e)
    reaction_step(r,e,[p["participant_id"] for p in parts],("100 °C",Q(100,"°C")),("2 h",Q(2,"h")))
    r["special_conditions"]["reactor"]=F("sealed tube","sealed tube",e)
    outcome(r,"yield_3a","yield",[pmid],e,"0%",Q(0,"%"),basis=None,yield_type=None)
    outcome(r,"control_observation","other",[],e,result_note,result_note)
    r["review"]["unresolved"]=["The measurement basis for the 0% desired-product result is not reported.","Conversion, recovery, and identities/amounts of any mixture or cycloisomerization product are unreported."]
    set_classification(r,"not_applicable" if extra else "unresolved",e); add_r(r)

# Scheme 5 is a proposal, not a sequence of executed experiments.
e=ev(4); r=base_reaction("scheme5_proposed_mechanism",e,"proposed_mechanism",["other"],"Scheme 5","Mechanistic proposal for the copper-catalyzed preparation of 4-(arylsulfonyl)-2,5-dihydrofurans 3","proposed radical copper catalytic cycle")
r["participants"]=[participant("allenol",f"{DOC}:mol_generic_allenol","reactant",e),participant("sulfinate",f"{DOC}:mol_generic_sulfinate","reactant",e),participant("product",f"{DOC}:mol_generic_product","product",e,product_role="proposed cycle product"),participant("copper",f"{DOC}:mol_CuOAc2","catalyst",e),participant("silver",f"{DOC}:mol_AgNO3","reagent",e)]
r["procedure"]["original_text"]=F("Author-proposed sequence: copper alkoxide/complex I, 5-endo-trig oxycupration to II, sulfinate oxidation to a sulfur radical, radical addition to III, reductive elimination, and Cu(I) reoxidation.","Scheme 5 proposed sequence",e)
r["review"]["unresolved"]=["Intermediates I-III are proposed and were not expanded into executed experiments; exact structures/atom mapping remain unresolved."]
set_classification(r,"cascade",e); add_r(r)

# Source-context records are intentionally not synthesis experiments.
bundle["source_context"] = [
    {"context_id":f"{DOC}:context_page1_metadata","document_id":DOC,"kind":"bibliographic_and_background","evidence_ids":ev(1),"note":"Title, authors, DOI, abstract, introduction, and online-SI notice; no current-paper concrete experiment arrow on page 1."},
    {"context_id":f"{DOC}:context_page4_characterization","document_id":DOC,"kind":"characterization","evidence_ids":ev(4),"note":"NMR, IR, and HRMS data for 3b are characterization linked to the 3b synthesis, not separate reactions."},
    {"context_id":f"{DOC}:context_page5_references","document_id":DOC,"kind":"references","evidence_ids":ev(5),"note":"References only; no current synthesis experiment, mass-spectrum table, molecular-weight table, or fragment-abundance block."},
    {"context_id":f"{DOC}:context_negative_spectra_inventory","document_id":DOC,"kind":"negative_inventory","evidence_ids":ev(1)+ev(2)+ev(3)+ev(4)+ev(5),"note":"No standalone MS peak list, molecular-weight table, or fragment-relative-abundance table was converted into a reaction."},
]

bundle["molecules"] = list(M.values())
bundle["reactions"] = R
assert len(R) == 40, len(R)
assert len({r["reaction_id"] for r in R}) == len(R)
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

# Raw visual inventory and stereo diagnostics are separate from formal records.
raw = {
    "document_id":DOC,
    "pages_inspected":[1,2,3,4,5],
    "record_count":len(R),
    "concrete_experiments":34,
    "inventory":{
        "page1":["title/metadata/abstract/introduction; no concrete current experiment"],
        "page2":["Scheme 1a,b,c,e cited prior work","Scheme 1d current generic scheme","Table 1 entries 1-13"],
        "page3":["Scheme 2 products 3a-i","Scheme 3a products 3j-p and failures 1m/1n","Scheme 3b furan 4a","Scheme 4 three controls"],
        "page4":["Scheme 5 proposed mechanism","Typical Procedure for 3b","3b characterization"],
        "page5":["references only"],
    },
    "deduplication":["Table 1 entry 2 and Scheme 2 product 3a merged","Scheme 2 product 3b and page-4 Typical Procedure merged"],
    "limitations":["No SI file supplied or associated","Compact drawings are candidate graph transcriptions, not verified structures","Pages beyond 1-5 were not claimed covered"],
}
(RUN/"raw_model/page_extraction.json").write_text(json.dumps(raw,ensure_ascii=False,indent=2)+"\n")
stereo = {}
for mid,m in M.items():
    smi=st.value(m["smiles"])
    if smi: stereo[mid]=st.stereo_inspect(smi)
(RUN/"raw_model/stereo_checks.json").write_text(json.dumps(stereo,ensure_ascii=False,indent=2)+"\n")

coverage={"document_id":DOC,"semantic_accuracy_verified":False,"pages":[
 {"page":1,"state":"attempted","regions":[{"region":"full page: title through introduction","finding":"Bibliographic metadata, abstract, introduction, and SI-availability notice; no current concrete experiment arrow.","evidence_ids":ev(1)}],"record_ids":[],"unresolved":["Advertised online SI was not supplied and is not associated."]},
 {"page":2,"state":"attempted","regions":[
   {"region":"upper-left Scheme 1","finding":"Four cited prior-work arrows (1a,b,c,e) and one current generic transformation (1d).","evidence_ids":ev(2)},
   {"region":"upper-right prose","finding":"Optimization discussion identifies 1a/2a, 100 °C, MeCN, sealed tube, 10 mol% Cu(OAc)2 and 200 mol% AgNO3; 3a 68%.","evidence_ids":ev(2)},
   {"region":"lower Table 1 including footnotes","finding":"All 13 optimization rows inventoried separately; entry 2 is the same 3a experiment repeated in Scheme 2.","evidence_ids":ev(2)}],
  "record_ids":[r["reaction_id"] for r in R if any(x in r["reaction_id"] for x in ("scheme1","table1"))],"unresolved":["Complete workups and atmospheres are absent from Table 1."]},
 {"page":3,"state":"attempted","regions":[
   {"region":"upper-left Scheme 2","finding":"Products 3a-i, isolated yields, and times; 3a merged with Table 1 entry 2.","evidence_ids":ev(3)},
   {"region":"upper/right Scheme 3a-b","finding":"Products 3j-p; nonproductive 1m/1n; prolonged-time aromatic furan 4a.","evidence_ids":ev(3)},
   {"region":"lower-right Scheme 4","finding":"Three mechanistic controls: no silver, no copper, and TEMPO inhibition.","evidence_ids":ev(3)},
   {"region":"lower-left and lower-right prose","finding":"Narrative limitations, radical interpretation, and start of proposed mechanism.","evidence_ids":ev(3)}],
  "record_ids":[r["reaction_id"] for r in R if any(x in r["reaction_id"] for x in ("scope_","nonproductive","furan_4a","control_","table1_entry2"))],"unresolved":["No structures/amounts are reported for messy products or detected cycloisomerization; no product was invented."]},
 {"page":4,"state":"attempted","regions":[
   {"region":"upper-left Scheme 5 and adjacent mechanism prose","finding":"Author-proposed catalytic cycle; not treated as executed intermediates.","evidence_ids":ev(4)},
   {"region":"lower-left Experimental Section continuing upper-right","finding":"Complete supplied Typical Procedure for 3b through chromatography; associated only with the 3b scope record.","evidence_ids":ev(4)},
   {"region":"upper-right characterization block","finding":"3b NMR/IR/HRMS retained as source context, not a new synthesis.","evidence_ids":ev(4)}],
  "record_ids":[f"{DOC}:scheme5_proposed_mechanism",f"{DOC}:scheme_scope_3b"],"unresolved":["Procedure gives TLC endpoint while Scheme 2 supplies 3 h; both are retained without treating them as a conflict."]},
 {"page":5,"state":"attempted","regions":[{"region":"both reference columns","finding":"References [17]-[29] only; no current-paper reaction arrow, procedure, or characterization table.","evidence_ids":ev(5)}],"record_ids":[],"unresolved":[]}
]}
(RUN/"coverage/pages_0001_0005.json").write_text(json.dumps(coverage,ensure_ascii=False,indent=2)+"\n")
