#!/usr/bin/env python3
"""Build source-specific v3.1 records for doc_1ba437add780cf5e, PDF pages 1-6."""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_1ba437add780cf5e"
TITLE = "Diastereodivergent Synthesis of β-Amino Alcohols by Dual-Metal-Catalyzed Coupling of Alkoxyallenes with Aldimine Esters"
DATASET = ROOT / "runs/current/dataset.json"


def T(page): return f"{DOC}:page{page}:text"
def I(page): return f"{DOC}:page{page}:image"
def ev(*pages):
    return [x for p in pages for x in (I(p), T(p))]


def F(raw=None, normalized=None, evidence=(), method="original_extraction", status=None,
      derivation=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def derived(normalized, evidence, derivation, raw="source structure/reaction drawing"):
    return F(raw, normalized, evidence, method="structure_derivation", derivation=derivation)


def q(number, unit, qualifier=None, raw_text=None):
    x = {"number": number, "unit": unit}
    if qualifier: x["qualifier"] = qualifier
    if raw_text: x["raw_text"] = raw_text
    return x


dataset = json.loads(DATASET.read_text())
doc = dataset["documents"][0]
assert doc["document_id"] == DOC
doc["title"] = TITLE
# Regeneration cleanup for the first machine draft: the scale-up produces the
# same source-identified 3ha molecule as Table 2, so remove that obsolete
# duplicate registry object before standard merge overwrites all reaction IDs.
if "standard" in dataset:
    dataset["standard"]["molecules"] = [m for m in dataset["standard"]["molecules"] if m["molecule_id"] != "m_3ha_scaled"]
DATASET.write_text(json.dumps(dataset, ensure_ascii=False, indent=2) + "\n")

bundle = st.empty(dataset["documents"])
molecules = {}
reactions = []


def mol(mid, label, pages, name=None, smiles=None, chirality="unresolved",
        geometry="unresolved", representation=None, assignment=None, axial=None,
        mixture=None, note=None):
    evidence = ev(*pages)
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(name or label, name or label, evidence)
    m["names"]["english"] = F(name or label, name or label, evidence)
    m["cas_numbers"] = F()
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smiles, evidence, "Candidate SMILES transcribed conservatively from the supplied source drawing; parsing does not verify identity.")
        diag = st.stereo_inspect(smiles)
        rep = representation or ("partial" if any(p["specified"] == "Unspecified" for p in diag["potential_stereo"]) else "encoded")
        m["stereochemistry"]["identity_key"] = derived(diag["stereo_graph_key"], evidence, "Computed from the candidate isomeric molecular graph.")
        if rep == "encoded":
            m["stereochemistry"]["isomeric_smiles"] = derived(diag["canonical_isomeric_smiles"], evidence, "Canonicalized from the candidate source-derived graph.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F()
        rep = representation or "not_encoded"
    stereoder = "Status derived from the explicit source label/drawing; missing stereomarks are not interpreted as racemic or achiral."
    m["stereochemistry"]["chirality_status"] = F("source label/drawing", chirality, evidence, method="inference", derivation=stereoder)
    m["stereochemistry"]["geometry_status"] = F("source label/drawing", geometry, evidence, method="inference", derivation=stereoder)
    m["stereochemistry"]["representation_status"] = F("candidate encoding", rep, evidence, method="inference", derivation="Compares the source's stereochemical content with the candidate representation.")
    if assignment is None and chirality in ("specified", "partially_specified"):
        assignment = "source-drawn or source-named stereochemical element; exact full assignment is not completely encoded"
    if assignment is not None:
        m["stereochemistry"]["assignment"] = F(assignment, assignment, evidence)
    if axial is not None:
        m["stereochemistry"]["axial"] = F(axial, axial, evidence)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = F(mixture, mixture, evidence)
    if note:
        m["identity_note"] = F(note, note, evidence)
    molecules[mid] = m
    return mid


def ensure(mid, label=None, pages=(2,), **kw):
    if mid not in molecules:
        mol(mid, label or mid, pages, **kw)
    return mid


def generic_mol(mid, label, pages, fixed, variables, reason):
    evidence = ev(*pages)
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(label, label, evidence)
    m["names"]["english"] = F(label, label, evidence)
    m["representation_type"] = "generic_template"
    m["smiles"] = F()
    m["generic_template"] = {
        "fixed_structure": F(fixed, fixed, evidence),
        "variables": [{
            "label": v,
            "attachment_points": F("shown attachment in source drawing", "shown attachment in source drawing", evidence),
            "definition": F(desc, desc, evidence),
            "allowed_values": F(),
            "constraints": F(),
        } for v, desc in variables],
        "instance_mapping": [],
        "reason": F(reason, reason, evidence),
    }
    for key, value in (("chirality_status", "unresolved"), ("geometry_status", "unresolved"), ("representation_status", "not_encoded")):
        m["stereochemistry"][key] = F("generic R-group drawing", value, evidence, method="inference", derivation=reason)
    molecules[mid] = m
    return mid


def fact_charge(raw, normalized, evidence): return F(raw, normalized, evidence)


def addp(r, mid, role, evidence, pid=None, stage="reaction", amount=None, equiv=None,
         molpct=None, volume=None, mass=None, ref=None, product_role=None):
    pid = pid or "p_" + re.sub(r"[^A-Za-z0-9]+", "_", mid).strip("_")
    p = st.new_participant(pid, mid, role, evidence, stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, evidence)
    for key, val in (("amount", amount), ("equivalents", equiv), ("mol_percent", molpct), ("volume", volume), ("mass", mass)):
        if val is not None:
            raw = val.get("raw_text") or f"{val.get('number')} {val.get('unit')}"
            p["charges"][key] = fact_charge(raw, val, evidence)
    if equiv is not None or molpct is not None:
        p["equivalence_reference"] = F(ref, ref, evidence)
    r["participants"].append(p)
    return pid


def add_step(r, sid, evidence, pids, operation, temperature=None, time=None, stage="reaction",
             atmosphere=None, pressure=None):
    s = st.new_step(sid, evidence, stage)
    s["participant_ids"] = list(pids)
    s["operation"] = F(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = F(temperature if isinstance(temperature, str) else temperature.get("raw_text", str(temperature)), temperature, evidence)
    if time is not None:
        s["time"] = F(time if isinstance(time, str) else time.get("raw_text", str(time)), time, evidence)
    if atmosphere is not None: s["atmosphere"] = F(atmosphere, atmosphere, evidence)
    if pressure is not None: s["pressure"] = F(pressure, pressure, evidence)
    r["steps"].append(s)


def outcome(r, oid, kind, targets, evidence, raw=None, normalized=None, yield_type=None,
            basis=None, scope=None, value_fact=None, basis_fact=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = value_fact or F(raw, normalized, evidence)
    if yield_type is not None: o["yield_type"] = F(yield_type, yield_type, evidence)
    if basis is not None: o["measurement_basis"] = F(basis, basis, evidence)
    if basis_fact is not None: o["measurement_basis"] = basis_fact
    if scope is not None: o["scope"] = F(scope, scope, evidence)
    if kind in ("yield", "product_amount"): o["product_role"] = F("product", "product", evidence)
    r["outcomes"].append(o)


TAX = {"C-H": 1, "C-C": 2, "C-O": 3, "C-Br": 9, "DB-C-C": 13,
       "O-H": 24, "N-H": 29, "C-N": 5}


def classify_partial(r, codes, evidence, limitation):
    changes, labels = [], []
    for i, code in enumerate(codes, 1):
        prefix = code.startswith("DB-")
        a, b = code.removeprefix("DB-").split("-")
        cid = f"b{i}"
        ch = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
              "before_order": 1 if prefix else 0, "after_order": 2 if prefix else 1,
              "change_type": "order_increased" if prefix else ("hydrogen_gain" if b == "H" else "formed")}
        if b == "H": ch["hydrogen_count"] = 2 if code == "N-H" else 1
        changes.append(ch)
        labels.append({"category_number": TAX[code], "code": code, "supporting_change_ids": [cid]})
    r["classification"]["bond_changes"] = [derived(c, evidence, limitation) for c in changes]
    r["classification"]["labels"] = derived(labels, evidence, limitation)
    r["classification"]["coverage"] = F("partial source-level transformation", "partial", evidence, method="inference", derivation=limitation)
    r["classification"]["limitations"] = [limitation, "No complete atom mapping is asserted; coproducts and additional hydrogen/bond-order changes may be absent."]


def classify_na(r, evidence, reason):
    r["classification"]["labels"] = F(None, [], evidence, method="inference", derivation=reason)
    r["classification"]["coverage"] = F(reason, "not_applicable", evidence, method="inference", derivation=reason)
    r["classification"]["limitations"] = [reason]


def rxn(suffix, pages, record_type="specific_experiment", categories=None, name=None,
        transform=None, location=None, general=True, unresolved=None):
    evidence = ev(*pages)
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories, categories, evidence) if categories is not None else F()
    r["location"] = F(location or f"PDF page(s) {', '.join(map(str, pages))}", {"pdf_pages": list(pages), "region": location}, evidence)
    if name: r["reaction_name"] = F(name, name, evidence)
    if transform: r["transformation_type"] = F(transform, transform, evidence)
    if general and record_type == "specific_experiment": r["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": unresolved or ["Supporting Information was not supplied; complete procedure and structure confirmation remain unavailable."]}
    r["completeness"] = "partial"
    reactions.append(r)
    return r


# Core molecules, generic families, catalysts, solvents, and workup agents.
generic_mol("m_generic_prior", "prior-work carbon nucleophile family", (1,),
            "azlactone/oxindole/pyrazolone + alkoxyallene generic structures",
            [("R", "source-defined substituent")], "Several cited substrates are grouped under one prior-work arrow.")
generic_mol("m_generic_product", "generic β-amino alcohol products", (1,),
            "branched allylic β-amino alcohol framework", [("R1", "aldimine substituent"), ("R2", "ester substituent"), ("R3", "alkoxyallene substituent")],
            "Scheme 1b is an R-group family, not a concrete molecule.")
generic_mol("m_generic_allene", "alkoxyallene family 1", (1, 4, 5), "R1O-CH=C=CH2",
            [("R1", "alkoxy substituent")], "R-group substrate family in the source.")
generic_mol("m_generic_aldimine", "aldimine ester family 2", (1, 4, 5), "Ar-CH=N-C(R2)-CO2R3",
            [("R2", "alpha substituent"), ("R3", "ester substituent")], "R-group substrate family in the source.")

mol("m_1a", "1a", (2, 3, 4, 5, 6), "benzyloxyallene 1a", "C=C=COCc1ccccc1",
    chirality="achiral", geometry="not_applicable", representation="partial",
    axial=[{"axis": "allene C=C=C", "assignment": "not stereogenic because one terminus is CH2"}],
    note="Candidate graph has a terminal CH2 allene; the local allene axis is not stereogenic on that structural basis.")
mol("m_2a", "2a", (2, 3, 4, 5), "p-fluorobenzylidene alanine methyl ester 2a",
    "COC(=O)C(C)N=Cc1ccc(F)cc1", chirality="unspecified", geometry="unspecified", representation="partial",
    note="Candidate graph retains an unassigned alpha stereocenter and unassigned imine geometry.")

for mid, label, name, smi in [
    ("m_Et3N", "Et3N", "triethylamine", "CCN(CC)CC"), ("m_THF", "THF", "tetrahydrofuran", "C1CCOC1"),
    ("m_PhMe", "PhMe", "toluene", "Cc1ccccc1"), ("m_DCM", "DCM", "dichloromethane", "ClCCl"),
    ("m_dioxane", "dioxane", "1,4-dioxane", "C1COCCO1"), ("m_CPME", "CPME", "cyclopentyl methyl ether", "COC1CCCC1"),
    ("m_citric", "aqueous citric acid (10%)", "citric acid", "OC(=O)CC(O)(CC(=O)O)C(=O)O"),
    ("m_AgOAc", "AgOAc", "silver acetate", "CC(=O)[O-].[Ag+]"),
    ("m_CuSalt", "Cu(CH3CN)4PF6", "tetrakis(acetonitrile)copper(I) hexafluorophosphate", None),
]:
    mol(mid, label, (2, 4, 5), name, smi, chirality="achiral" if smi else "unresolved",
        geometry="not_applicable" if smi else "unresolved", representation="encoded" if smi else "not_encoded")

for lig in ["L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11", "L12", "L13", "ent-L1", "ent-L3", "ent-L13"]:
    assignment = "source-specified chiral ligand identity" if lig not in ("L9", "L10") else "source-specified ligand identity"
    mol(f"m_{lig.replace('-', '_')}", lig, (2, 3, 4, 5, 6), lig, chirality="specified", geometry="not_applicable",
        representation="not_encoded", assignment=assignment)
for lig in ["L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11", "L12", "L13", "ent-L13", "ent-L3"]:
    mol(f"m_Pd_{lig.replace('-', '_')}", f"{lig}-Pd(allyl)BF4", (2, 3, 4, 5, 6), f"{lig}-Pd(allyl)BF4",
        chirality="specified", geometry="not_applicable", representation="not_encoded", assignment="source-specified chiral Pd catalyst")

for mid, label, name, smi, page in [
    ("m_Boc2O", "(Boc)2O", "di-tert-butyl dicarbonate", "CC(C)(C)OC(=O)OC(=O)OC(C)(C)C", 6),
    ("m_LiAlH4", "LiAlH4", "lithium aluminium hydride", "[Li+].[AlH4-]", 6),
    ("m_O3", "O3", "ozone", "[O-][O+]=O", 6), ("m_NaBH4", "NaBH4", "sodium borohydride", "[Na+].[BH4-]", 6),
    ("m_MeOH", "MeOH", "methanol", "CO", 6), ("m_PdC", "Pd/C", "palladium on carbon", None, 6),
    ("m_H2", "H2", "hydrogen", "[H][H]", 6), ("m_NBS", "N-bromosuccinimide", "N-bromosuccinimide", "O=C1CCC(=O)N1Br", 6),
    ("m_Na2CO3", "Na2CO3", "sodium carbonate", "[Na+].[Na+].[O-]C([O-])=O", 6),
    ("m_MeCN", "MeCN", "acetonitrile", "CC#N", 6),
]:
    if mid == "m_LiAlH4":
        smi = None
    mol(mid, label, (page,), name, smi, chirality="achiral" if smi else "unresolved",
        geometry="not_applicable" if smi else "unresolved", representation="encoded" if smi else "not_encoded")

# Concrete product stereoisomers used across tables and schemes.
base_product_smiles = "COC(=O)C(N)(C)C(OCC1=CC=CC=C1)C=C"
mol("m_3aa", "(2S,3S)-3aa", (2, 3, 4, 6), "(2S,3S)-3aa", base_product_smiles,
    chirality="specified", geometry="not_applicable", representation="partial", assignment="(2S,3S); anti")
mol("m_4aa", "(2S,3R)-4aa", (2, 3, 5, 6), "(2S,3R)-4aa", base_product_smiles,
    chirality="specified", geometry="not_applicable", representation="partial", assignment="(2S,3R); syn")
mol("m_3aa_RR", "(2R,3R)-3aa", (3,), "(2R,3R)-3aa", base_product_smiles,
    chirality="specified", geometry="not_applicable", representation="partial", assignment="(2R,3R); anti")
mol("m_4aa_RS", "(2R,3S)-4aa", (3,), "(2R,3S)-4aa", base_product_smiles,
    chirality="specified", geometry="not_applicable", representation="partial", assignment="(2R,3S); syn")

# Scheme 1: one cited family and one current generic transformation.
r = rxn("scheme1_prior_work", (1,), "literature_example", None,
        "Pd-catalyzed diastereoselective Csp3-Csp3 coupling", "addition of carbon nucleophiles to alkoxyallenes",
        "Scheme 1a", general=False, unresolved=["Cited azlactone, oxindole, and pyrazolone examples are grouped by the source; conditions and individual outcomes are not supplied here."])
addp(r, "m_generic_prior", "reactant", ev(1)); addp(r, "m_generic_allene", "reactant", ev(1), pid="p_allene")
classify_partial(r, ["C-C"], ev(1), "Prior-work family classification is limited to the source's Csp3-Csp3 coupling description.")

r = rxn("scheme1_current_generic", (1,), "generic_scheme", None,
        "dual-metal-catalyzed diastereodivergent coupling", "alkoxyallene/aldimine-ester coupling", "Scheme 1b", general=False,
        unresolved=["R-group family only; no concrete experiment, yield, or unique structures asserted."])
addp(r, "m_generic_allene", "reactant", ev(1)); addp(r, "m_generic_aldimine", "reactant", ev(1), pid="p_aldimine")
addp(r, "m_generic_product", "product", ev(1), pid="p_products", product_role="syn/anti generic products")
classify_partial(r, ["C-C", "N-H"], ev(1), "Class-level changes inferred from the generic coupling drawing and acid-unmasked primary amine product.")


# Table 1 optimization, 22 independent rows.
table1 = [
 (1,"L4","L1","Cu","THF","2:1",57,"98","97",None),(2,"L5","L1","Cu","THF","1.4:1",63,"99","97",None),
 (3,"L6","L1","Cu","THF","3.7:1",66,">99","99",None),(4,"L7","L1","Cu","THF","1:2.2",54,"97","98",None),
 (5,"L8","L1","Cu","THF","1:1.6",70,"81","97",None),(6,"L9","L1","Cu","THF","1.1:1",69,"82","98",None),
 (7,"L10","L1","Cu","THF","1:2",65,">99","97",None),(8,"L11","L1","Cu","THF","1:4.7",65,">99",">99",None),
 (9,"L12","L1","Cu","THF","7:1",78,">99",">99",None),(10,"L13","L1","Cu","THF","14:1",80,">99",None,"conflict70_71"),
 (11,"ent-L13","L1","Cu","THF","1:3.9",64,">99",">99",None),(12,"L1","L1","Cu","THF","1:4.1",73,">99",">99",None),
 (13,"L2","L2","Cu","THF","1:4.5",64,"99",">99",None),(14,"L3","L3","Cu","THF","1:5.2",60,">99",">99",None),
 (15,"L3","L3","Ag","THF","1:12",42,None,">99",None),(16,"L3","L3","Ag","PhMe","1:14",72,None,">99",None),
 (17,"L3","L3","Ag","DCM","1:2.5",54,"41","90",None),(18,"L3","L3","Ag","dioxane","1:14",48,None,">99",None),
 (19,"L3","L3","Ag","CPME","1:18",80,None,">99",None),(20,"L3","L3","Ag","CPME","1:20",82,None,">99",73),
 (21,"L13",None,None,"THF",None,"<5",None,None,None),(22,"L3",None,None,"CPME",None,"<5",None,None,None),
]

for entry, la, lb, metal, solvent, ratio, yld, eeanti, eesyn, isolated in table1:
    pages = (2, 3) if entry in (10, 20) else (2,)
    evidence = ev(*pages)
    r = rxn(f"table1_entry{entry}", pages, categories=["optimization"],
            name="optimization of alkoxyallene/aldimine-ester coupling", transform="diastereodivergent coupling",
            location=f"Table 1, entry {entry}" + ("; Scheme 2 repeated outcome" if entry in (10,20) else ""),
            unresolved=["SI was not supplied; complete addition order and preparative workup are unavailable."] +
                       (["Table 1 gives 70% isolated for entry 10, whereas prose/Scheme 2 give 71%; both are retained without choosing."] if entry == 10 else []))
    p1 = addp(r, "m_1a", "reactant", evidence, amount=q(0.1,"mmol"), equiv=q(1,"equiv"), ref="m_1a")
    p2 = addp(r, "m_2a", "reactant", evidence, amount=q(0.2,"mmol"), equiv=q(2,"equiv"), ref="m_1a")
    addp(r, "m_3aa", "product", evidence, product_role="anti product")
    addp(r, "m_4aa", "product", evidence, product_role="syn product")
    pdmid = f"m_Pd_{la.replace('-', '_')}"; addp(r, pdmid, "catalyst", evidence, molpct=q(4,"mol%"), ref="m_1a")
    if lb:
        addp(r, f"m_{lb.replace('-', '_')}", "ligand", evidence, molpct=q(5 if entry==20 else 6,"mol%"), ref="m_1a")
    if metal == "Cu": addp(r, "m_CuSalt", "catalyst", evidence, molpct=q(5,"mol%"), ref="m_1a")
    if metal == "Ag": addp(r, "m_AgOAc", "catalyst", evidence, molpct=q(6 if entry==20 else 5,"mol%"), ref="m_1a")
    addp(r, "m_Et3N", "reagent", evidence, molpct=q(200,"mol%"), ref="m_1a")
    smid = f"m_{solvent}"; addp(r, smid, "solvent", evidence, volume=q(0.5,"mL"))
    r["reaction_scale"] = F("1a (0.1 mmol)", {"reference_molecule_id":"m_1a","amount":q(0.1,"mmol")}, evidence)
    r["solvent_conditions"]["composition"] = F(solvent, [smid], evidence)
    r["solvent_conditions"]["total_volume"] = F("0.5 mL", q(0.5,"mL"), evidence)
    add_step(r, "reaction", evidence, [p1,p2], "dual-metal-catalyzed coupling", q(30,"°C"), q(48,"h"))
    pc = addp(r, "m_citric", "workup_agent", evidence, stage="workup")
    add_step(r, "citric_acid_workup", evidence, [pc], "treat with aqueous citric acid (10%)", None, q(1,"h"), stage="workup")
    ynorm = q(5,"%", "lt", "<5%") if yld == "<5" else q(yld,"%")
    outcome(r, "crude_yield", "yield", ["m_3aa","m_4aa"], evidence, raw=f"{yld}%", normalized=ynorm,
            yield_type="combined/mixture yield", basis="1H NMR of crude products")
    if isolated == "conflict70_71":
        vf = F(["70% isolated (Table 1)", "71% isolated (prose/Scheme 2)"], None, evidence, status="conflict",
               alternatives=[F("70%",q(70,"%"),ev(2)),F("71%",q(71,"%"),ev(2,3))])
        outcome(r, "isolated_yield_major", "yield", ["m_3aa"], evidence, value_fact=vf, yield_type="isolated", basis="isolated")
    elif isolated is not None:
        outcome(r, "isolated_yield_major", "yield", ["m_4aa" if entry==20 else "m_3aa"], evidence,
                raw=f"{isolated}%", normalized=q(isolated,"%"), yield_type="isolated", basis="isolated")
    if ratio:
        outcome(r, "anti_syn_ratio", "selectivity", ["m_3aa","m_4aa"], evidence, raw=ratio, normalized=ratio,
                basis="1H NMR of crude products", scope={"metric":"dr","order":["anti 3aa","syn 4aa"]})
    for oid, target, value in (("ee_anti","m_3aa",eeanti),("ee_syn","m_4aa",eesyn)):
        if value:
            norm=q(float(value.replace(">","")),"%", "gt" if value.startswith(">") else None, value+"%")
            outcome(r, oid, "selectivity", [target], evidence, raw=value+"% ee", normalized=norm,
                    basis="chiral HPLC", scope={"metric":"ee"})
    classify_partial(r, ["C-C","N-H"], evidence, "Partial net classification inferred from the coupling product and acid hydrolysis of the aldimine auxiliary.")


# The two catalyst-inverted Scheme 2 experiments not present as Table 1 rows.
for suffix, product, catalysts, yld, dr in [
    ("scheme2_2R3R_3aa","m_3aa_RR",("m_Pd_ent_L13","m_ent_L1","m_CuSalt"),72,"17:1"),
    ("scheme2_2R3S_4aa","m_4aa_RS",("m_Pd_ent_L3","m_ent_L3","m_AgOAc"),76,"15:1"),
]:
    r = rxn(suffix,(3,),categories=["substrate_scope"],name="stereodivergent coupling",transform="diastereodivergent coupling",location="Scheme 2",
            unresolved=["Scheme 2 does not state the yield or dr measurement basis; complete catalyst charges and procedure are in unavailable SI."])
    addp(r,"m_1a","reactant",ev(3)); addp(r,"m_2a","reactant",ev(3),pid="p_2a"); addp(r,product,"product",ev(3),product_role="major stereoisomer")
    for j,mid in enumerate(catalysts): addp(r,mid,"catalyst" if "Pd_" in mid or mid.endswith("Salt") or mid=="m_AgOAc" else "ligand",ev(3),pid=f"p_cat{j}")
    add_step(r,"reaction",ev(3),["p_m_1a","p_2a"],"couple 1a and 2a with the displayed catalyst pair")
    outcome(r,"yield","yield",[product],ev(3),raw=f"{yld}%",normalized=q(yld,"%"))
    outcome(r,"dr","selectivity",[product],ev(3),raw=dr,normalized=dr,scope={"metric":"dr"})
    outcome(r,"ee","selectivity",[product],ev(3),raw=">99% ee",normalized=q(99,"%","gt",">99%"),scope={"metric":"ee"})
    classify_partial(r,["C-C","N-H"],ev(3),"Partial net classification inferred from Scheme 2 product connectivity and primary amine formation.")


def ensure_substrate(label, page, ch="unspecified", assignment=None):
    mid="m_"+label.replace("-","_").replace("(","").replace(")","")
    return ensure(mid,label,(page,),name=label,chirality=ch,geometry="unspecified",representation="not_encoded",assignment=assignment)


def ensure_product(label, page, assignment, ch="specified"):
    mid="m_"+label.replace("-","_").replace("(","").replace(")","").replace(",","_")
    return ensure(mid,label,(page,),name=label,chirality=ch,geometry="not_applicable",representation="not_encoded",assignment=assignment)


# Table 2 anti-selective scope. Data are read row-by-row from the source image.
anti = [
 ("3ab",66,"16:1",">99"),("3ac",65,"13:1",">99"),("3ad",54,"13:1",">99"),("3ae",69,"10:1",">99"),("3af",77,"16:1","96"),
 ("3ag",70,"20:1",">99"),("3ah",62,"15:1",">99"),("3ai",78,"14:1",">99"),("3aj",70,">20:1",">99"),("3ak",48,"19:1",">99"),
 ("3al",75,">20:1",None),("3am",65,"14:1",None),("3an",62,"12:1",None),("3ao",69,"12:1",">99"),("3ap",51,">20:1",None),
 ("3aa",71,"14:1",">99"),("3ba",71,"19:1",">99"),("3ca",43,">20:1",">99"),("3da",62,"14:1",">99"),("3ea",56,"9:1",">99"),
 ("3fa",85,"16:1",">99"),("3ga",72,">20:1",">99"),("3ha",70,"11:1","99"),("3ia",62,"16:1",">99"),("3ja",None,"4:1",">99"),
]

for label,yld,dr,ee in anti:
    if label == "3aa":
        product="m_3aa"; s1="m_1a"; s2="m_2a"
    elif label.startswith("3a"):
        letter=label[-1]; s1=ensure_substrate("1"+letter,4,"partially_specified" if letter in "lmnop" else "unspecified"); s2="m_2a"
        product=ensure_product(f"(2S,3S)-{label}",4,"(2S,3S); anti" + (" plus source-drawn substrate stereochemistry" if letter in "lmnop" else ""),"partially_specified" if letter in "lmnop" else "specified")
    else:
        letter=label[1]; s1="m_1a"; s2=ensure_substrate("2"+letter,4,"unspecified")
        product=ensure_product(f"(2S,3S)-{label}",4,"(2S,3S); anti")
    evidence=ev(3,4) if label in ("3ha","3ia","3ja") else ev(4)
    unresolved=["Substrate/product structures are retained from the compact drawing without a complete verified SMILES; SI was not supplied."]
    if label=="3ja": unresolved.append("Table 2 shows 83% isolated yield, while page-3 prose states 56%; both are retained as a conflict.")
    r=rxn(f"table2_{label}",(3,4) if label in ("3ha","3ia","3ja") else (4,),categories=["substrate_scope"],
          name="anti-selective coupling",transform="alkoxyallene/aldimine-ester coupling",location=f"Table 2, product {label}",unresolved=unresolved)
    addp(r,s1,"reactant",evidence,amount=q(0.1,"mmol"),equiv=q(1,"equiv"),ref=s1)
    addp(r,s2,"reactant",evidence,pid="p_substrate2",amount=q(0.2,"mmol"),equiv=q(2,"equiv"),ref=s1)
    addp(r,product,"product",evidence,product_role="major anti product")
    addp(r,"m_Pd_L13","catalyst",evidence,molpct=q(4,"mol%"),ref=s1)
    addp(r,"m_L1","ligand",evidence,molpct=q(6,"mol%"),ref=s1)
    addp(r,"m_CuSalt","catalyst",evidence,molpct=q(5,"mol%"),ref=s1)
    addp(r,"m_Et3N","reagent",evidence,molpct=q(200,"mol%"),ref=s1)
    addp(r,"m_THF","solvent",evidence,volume=q(0.5,"mL"))
    r["reaction_scale"]=F("alkoxyallene 1 (0.1 mmol)",{"reference_molecule_id":s1,"amount":q(0.1,"mmol")},evidence)
    r["solvent_conditions"]["composition"]=F("THF",["m_THF"],evidence); r["solvent_conditions"]["total_volume"]=F("0.5 mL",q(0.5,"mL"),evidence)
    add_step(r,"reaction",evidence,[r["participants"][0]["participant_id"],"p_substrate2"],"anti-selective dual-metal coupling",q(30,"°C"),q(48,"h"))
    pc=addp(r,"m_citric","workup_agent",evidence,stage="workup"); add_step(r,"workup",evidence,[pc],"treat with aqueous citric acid (10%)",None,q(1,"h"),stage="workup")
    if label=="3ja":
        vf=F(["83% (Table 2)","56% (page-3 prose)"],None,evidence,status="conflict",alternatives=[F("83%",q(83,"%"),ev(4)),F("56%",q(56,"%"),ev(3))])
        bf=F(["isolated yields are given (Table 2)","prose basis unreported"],None,evidence,status="conflict",alternatives=[F("isolated","isolated",ev(4)),F()])
        outcome(r,"yield","yield",[product],evidence,value_fact=vf,yield_type="source-conflicted",basis_fact=bf)
    else:
        outcome(r,"yield","yield",[product],evidence,raw=f"{yld}%",normalized=q(yld,"%"),yield_type="isolated",basis="isolated")
    outcome(r,"dr","selectivity",[product],evidence,raw=dr,normalized=dr,basis="1H NMR of crude products",scope={"metric":"dr","major":"anti"})
    if ee:
        outcome(r,"ee","selectivity",[product],evidence,raw=ee+"% ee",normalized=q(float(ee.replace(">","")),"%","gt" if ee.startswith(">") else None,ee+"%"),basis="chiral HPLC",scope={"metric":"ee"})
    classify_partial(r,["C-C","N-H"],evidence,"Partial net classification inferred from Table 2 product drawings and acid hydrolysis; exact atom mapping and coproducts are unresolved.")


# Table 3 syn-selective scope.
syn = [
 ("4ab",71,"20:1","99"),("4ac",72,"12:1","99"),("4ad",52,"16:1","99"),("4ae",66,"13:1",">99"),("4af",62,"14:1",">99"),
 ("4ag",64,"13:1",">99"),("4ah",41,"20:1",">99"),("4ai",62,"6:1","96"),("4aj",62,"8:1",None),("4ak",54,"4:1",None),
 ("4al",69,"18:1",None),("4am",63,">20:1",None),("4an",57,"17:1",None),
 ("4aa",73,"20:1",">99"),("4ba",63,"13:1",">99"),("4ca",40,"20:1","98"),("4da",58,"13:1",">99"),("4ea",50,"9:1",">99"),
 ("4fa",81,"17:1",">99"),("4ga",71,"14:1",">99"),("4ha",64,"18:1","96"),("4ia",63,"7:1","80"),("4ja",85,"2:1","95"),
]
for label,yld,dr,ee in syn:
    pages=(3,5) if label in ("4am","4an","4ia","4ja") else (5,)
    evidence=ev(*pages)
    if label=="4aa": product="m_4aa"; s1="m_1a"; s2="m_2a"
    elif label.startswith("4a"):
        letter=label[-1]
        if letter=="m": s1=ensure_substrate("1m / (R)-1x",5,"specified","(R)-1-phenylethanol-derived substrate")
        elif letter=="n": s1=ensure_substrate("1n / (S)-1x",5,"specified","(S)-1-phenylethanol-derived substrate")
        else: s1=ensure_substrate("1"+letter,5,"partially_specified" if letter in "ijkl" else "unspecified")
        s2="m_2a"
        assn={"m":"(2S,3R,R); syn","n":"(2S,3R,S); syn"}.get(letter,"(2S,3R); syn" + (" plus source-drawn substrate stereochemistry" if letter in "ijkl" else ""))
        product=ensure_product(assn.split("; syn")[0]+f"-{label}",5,assn,"partially_specified" if letter in "ijkl" else "specified")
    else:
        letter=label[1]; s1="m_1a"; s2=ensure_substrate("2"+letter,5,"unspecified"); product=ensure_product(f"(2S,3R)-{label}",5,"(2S,3R); syn")
    r=rxn(f"table3_{label}",pages,categories=["substrate_scope"],name="syn-selective coupling",transform="alkoxyallene/aldimine-ester coupling",location=f"Table 3, product {label}" + ("; Scheme 3 repeated outcome" if label in ("4am","4an") else ""),
          unresolved=["Compact source structures are not asserted as complete verified SMILES; SI was not supplied."])
    addp(r,s1,"reactant",evidence,amount=q(0.1,"mmol"),equiv=q(1,"equiv"),ref=s1)
    addp(r,s2,"reactant",evidence,pid="p_substrate2",amount=q(0.2,"mmol"),equiv=q(2,"equiv"),ref=s1)
    addp(r,product,"product",evidence,product_role="major syn product")
    addp(r,"m_Pd_L3","catalyst",evidence,molpct=q(4,"mol%"),ref=s1); addp(r,"m_L3","ligand",evidence,molpct=q(5,"mol%"),ref=s1)
    addp(r,"m_AgOAc","catalyst",evidence,molpct=q(6,"mol%"),ref=s1); addp(r,"m_Et3N","reagent",evidence,molpct=q(200,"mol%"),ref=s1)
    addp(r,"m_CPME","solvent",evidence,volume=q(0.5,"mL"))
    r["reaction_scale"]=F("alkoxyallene 1 (0.1 mmol)",{"reference_molecule_id":s1,"amount":q(0.1,"mmol")},evidence)
    r["solvent_conditions"]["composition"]=F("CPME",["m_CPME"],evidence); r["solvent_conditions"]["total_volume"]=F("0.5 mL",q(0.5,"mL"),evidence)
    add_step(r,"reaction",evidence,[r["participants"][0]["participant_id"],"p_substrate2"],"syn-selective dual-metal coupling",q(30,"°C"),q(48,"h"))
    pc=addp(r,"m_citric","workup_agent",evidence,stage="workup"); add_step(r,"workup",evidence,[pc],"treat with aqueous citric acid (10%)",None,q(1,"h"),stage="workup")
    outcome(r,"yield","yield",[product],evidence,raw=f"{yld}%",normalized=q(yld,"%"),yield_type="isolated",basis="isolated")
    outcome(r,"dr","selectivity",[product],evidence,raw=dr,normalized=dr,basis="1H NMR of crude products",scope={"metric":"dr","major":"syn"})
    if ee: outcome(r,"ee","selectivity",[product],evidence,raw=ee+"% ee",normalized=q(float(ee.replace(">","")),"%","gt" if ee.startswith(">") else None,ee+"%"),basis="chiral HPLC",scope={"metric":"ee"})
    classify_partial(r,["C-C","N-H"],evidence,"Partial net classification inferred from Table 3 product drawings and acid hydrolysis; exact atom mapping and coproducts are unresolved.")


# Scheme 3: six additional catalyst-controlled stereoisomers; two standard Ag/Pd rows are already Table 3 4am/4an.
scheme3 = [
 ("scheme3_4an_2S3SS","1n / (S)-1x","(2S,3S,S)-4an",64,"6:1",("m_Pd_L13","m_L1","m_CuSalt")),
 ("scheme3_4an_2R3RS","1n / (S)-1x","(2R,3R,S)-4an",67,">20:1",("m_Pd_ent_L13","m_ent_L1","m_CuSalt")),
 ("scheme3_4an_2R3SS","1n / (S)-1x","(2R,3S,S)-4an",58,">20:1",("m_Pd_ent_L3","m_ent_L3","m_AgOAc")),
 ("scheme3_4am_2S3SR","1m / (R)-1x","(2S,3S,R)-4am",78,">20:1",("m_Pd_L13","m_L1","m_CuSalt")),
 ("scheme3_4am_2R3RR","1m / (R)-1x","(2R,3R,R)-4am",65,"5:1",("m_Pd_ent_L13","m_ent_L1","m_CuSalt")),
 ("scheme3_4am_2R3SR","1m / (R)-1x","(2R,3S,R)-4am",55,">20:1",("m_Pd_ent_L3","m_ent_L3","m_AgOAc")),
]
for suffix,substrate_label,prod_label,yld,dr,cats in scheme3:
    s1=ensure_substrate(substrate_label,5,"specified","(R)/(S)-1-phenylethanol-derived stereocenter as printed")
    product=ensure_product(prod_label,5,prod_label.split("-")[0])
    r=rxn(suffix,(5,),categories=["substrate_scope"],name="stereodivergent synthesis with three stereocenters",transform="diastereodivergent coupling",location="Scheme 3",
          unresolved=["Scheme 3 does not report yield/dr measurement methods or detailed catalyst charges."])
    addp(r,s1,"reactant",ev(5)); addp(r,"m_2a","reactant",ev(5),pid="p_2a"); addp(r,product,"product",ev(5),product_role="displayed stereoisomer")
    for j,mid in enumerate(cats): addp(r,mid,"catalyst" if "Pd_" in mid or mid in ("m_CuSalt","m_AgOAc") else "ligand",ev(5),pid=f"p_cat{j}")
    add_step(r,"reaction",ev(5),[r["participants"][0]["participant_id"],"p_2a"],"couple chiral alkoxyallene with 2a using the displayed catalyst pair")
    outcome(r,"yield","yield",[product],ev(5),raw=f"{yld}%",normalized=q(yld,"%")); outcome(r,"dr","selectivity",[product],ev(5),raw=dr,normalized=dr,scope={"metric":"dr"})
    classify_partial(r,["C-C","N-H"],ev(5),"Partial net classification inferred from Scheme 3 coupling products; exact mapping is unresolved.")


# Scheme 4 is author-proposed stereochemical rationale, not an executed experiment.
r=rxn("scheme4_proposed_mechanism",(4,5,6),"proposed_mechanism",None,"rationalization of stereoselectivities","proposed dual-metal catalytic mechanism","Scheme 4 and adjacent discussion",general=False,
      unresolved=["Intermediates A-F are proposed catalyst-bound species, not isolated experiments; transition-state structures are not established."])
for mid in ("m_1a","m_2a","m_3aa","m_4aa","m_L1","m_L3","m_CuSalt","m_AgOAc","m_Pd_L13","m_Pd_L3"):
    addp(r,mid,"product" if mid in ("m_3aa","m_4aa") else "reactant" if mid in ("m_1a","m_2a") else "catalyst",ev(4,5,6),pid="p_"+mid)
r["classification"]["coverage"]=F("proposed mechanism; no executed net experiment", "unresolved", ev(4,5,6), method="inference", derivation="Mechanistic arrows are hypotheses rather than measured transformations.")
r["classification"]["limitations"]=["No executed experiment or complete atom mapping is asserted for proposed intermediates."]


# Scheme 5 transformations. The clear 5a drawing supports a conservative
# connectivity-only candidate; its specified stereochemistry remains partial.
mol("m_5a", "5a", (6,), "Boc-protected (2S,3S)-3aa", "COC(=O)C(NC(=O)OC(C)(C)C)(C)C(OCC1=CC=CC=C1)C=C",
    chirality="specified", geometry="not_applicable", representation="partial", assignment="(2S,3S); Boc-protected")
for label,assn in [("5a","(2S,3S); Boc-protected"),("5b","source-drawn major diastereomer"),("5c","source-drawn major diastereomer"),("5d","source-drawn major diastereomer"),("5e","source-drawn major diastereomer")]:
    ensure_product(label,6,assn,"partially_specified")

scheme5 = [
 ("5a","m_3aa",["m_Boc2O","m_Et3N","m_THF"],73,"(Boc)2O, Et3N, THF; 0 °C to rt; 15 h",["C-N"]),
 ("5b","m_5a",["m_LiAlH4","m_THF"],61,"LiAlH4, THF; 0 °C to rt; 12 h",["C-O","C-H","O-H"]),
 ("5c","m_5a",["m_O3","m_NaBH4","m_DCM","m_MeOH"],54,"O3 then NaBH4; DCM/MeOH 3:1; -78 °C to rt; 2 h",["C-O","C-H","O-H"]),
 ("5d","m_5a",["m_PdC","m_H2","m_MeOH"],78,"Pd/C (10 mol%), H2 (1 atm balloon), MeOH, rt, 5 h",["C-H","O-H"]),
 ("5e","m_5a",["m_NBS","m_Na2CO3","m_MeCN"],72,"NBS, Na2CO3, MeCN, rt, 15 h",["C-N","C-Br"]),
]
for label,start,reags,yld,cond,codes in scheme5:
    product="m_"+label; r=rxn(f"scheme5_{label}",(6,),categories=["derivatization"],name="synthetic transformation of coupling product",transform=cond.split(";")[0],location=f"Scheme 5, arrow {label}",
        unresolved=["Scheme 5 does not state yield/dr measurement methods, reagent amounts (except Pd/C loading/pressure), workup, or purification."])
    addp(r,start,"reactant",ev(6)); addp(r,product,"product",ev(6),product_role="major product")
    pids=[]
    for j,mid in enumerate(reags):
        kwargs={}
        if mid=="m_PdC": kwargs={"molpct":q(10,"mol%"),"ref":start}
        if mid=="m_H2": kwargs={}
        pids.append(addp(r,mid,"solvent" if mid in ("m_THF","m_DCM","m_MeOH","m_MeCN") else "catalyst" if mid=="m_PdC" else "reagent",ev(6),pid=f"p_r{j}",**kwargs))
    # Explicitly separate initial low-temperature state from warming/hold when the source gives a range.
    if label in ("5a","5b"):
        add_step(r,"low_temperature_addition",ev(6),[r["participants"][0]["participant_id"]]+pids,"combine at 0 °C",q(0,"°C"))
        add_step(r,"warm_and_hold",ev(6),[r["participants"][0]["participant_id"]]+pids,"warm to room temperature and hold","rt",q(15 if label=="5a" else 12,"h"))
    elif label=="5c":
        add_step(r,"ozonolysis_and_reduction",ev(6),[r["participants"][0]["participant_id"]]+pids,"ozonolysis then NaBH4 reduction; warm from -78 °C to room temperature",{"minimum":-78,"maximum":20,"unit":"°C","raw_text":"-78 °C to rt"},q(2,"h"))
        r["solvent_conditions"]["composition"]=F("3:1 dichloromethane/MeOH",["m_DCM","m_MeOH"],ev(6)); r["solvent_conditions"]["mixture_ratio"]=F("3:1",{"components":["m_DCM","m_MeOH"],"ratio":"3:1"},ev(6))
    elif label=="5d": add_step(r,"hydrogenation",ev(6),[r["participants"][0]["participant_id"]]+pids,"hydrogenate under a balloon","rt",q(5,"h"),pressure="1 atm balloon")
    else: add_step(r,"bromocyclization",ev(6),[r["participants"][0]["participant_id"]]+pids,"NBS-induced bromocyclization","rt",q(15,"h"))
    outcome(r,"yield","yield",[product],ev(6),raw=f"{yld}%",normalized=q(yld,"%")); outcome(r,"dr","selectivity",[product],ev(6),raw=">20:1",normalized=">20:1",scope={"metric":"dr"})
    if label=="5a": outcome(r,"ee","selectivity",[product],ev(6),raw=">99% ee",normalized=q(99,"%","gt",">99%"),scope={"metric":"ee"})
    classify_partial(r,codes,ev(6),f"Partial class-level changes inferred from the Scheme 5 {label} drawing; exact atoms and all order/cleavage changes are unresolved.")


# Page-6 prose experiments in the mycestericin syntheses. Scheme 6 itself is on page 7 and is outside claimed coverage.
for mid,label,name in [("m_2h","2h","oxazoline 2h"),
                       ("m_6a","6a","alcohol 6a"),("m_7","7","ketoalkene 7"),("m_8a","8a","metathesis product 8a"),
                       ("m_mycF","mycestericin F","mycestericin F"),("m_3hb","(2S,3R)-3hb","source-labelled syn-coupling product 3hb"),
                       ("m_mycG","mycestericin G","mycestericin G")]:
    mol(mid,label,(6,),name,chirality="specified" if any(x in label for x in ("2S","mycestericin")) else "unspecified",
        geometry="unspecified",representation="not_encoded",assignment=label if any(x in label for x in ("2S","mycestericin")) else None,
        note="Compact structure/prose identity retained; Scheme 6 is outside the claimed page range and no complete SMILES is asserted.")
mol("m_HGII","Hoveyda-Grubbs II catalyst",(6,),"Hoveyda-Grubbs II catalyst",chirality="unresolved",geometry="unresolved",representation="not_encoded")
# The scaled experiment reuses the same 3ha molecular identity as Table 2;
# extend its occurrence evidence instead of creating a scale-specific molecule.
for eid in ev(6):
    if eid not in molecules["m_2S_3S_3ha"]["occurrences"][0]["evidence_ids"]:
        molecules["m_2S_3S_3ha"]["occurrences"][0]["evidence_ids"].append(eid)

r=rxn("scheme6_scaled_3ha",(6,),categories=["scale_up"],name="scaled-up anti-selective coupling",transform="alkoxyallene/oxazoline coupling",location="page-6 prose introducing Scheme 6",
      unresolved=["The exact scale, charges, temperature, time, and yield measurement method are not stated on page 6; Scheme 6 is outside this session's covered pages."])
addp(r,"m_1a","reactant",ev(6)); addp(r,"m_2h","reactant",ev(6),pid="p_2h"); addp(r,"m_2S_3S_3ha","product",ev(6),product_role="anti product")
addp(r,"m_Pd_L13","catalyst",ev(6)); addp(r,"m_L1","ligand",ev(6),pid="p_L1"); addp(r,"m_CuSalt","catalyst",ev(6),pid="p_Cu")
add_step(r,"reaction",ev(6),["p_m_1a","p_2h"],"scaled-up coupling under the L13-Pd/L1-Cu system")
outcome(r,"yield","yield",["m_2S_3S_3ha"],ev(6),raw="60%",normalized=q(60,"%")); outcome(r,"ee","selectivity",["m_2S_3S_3ha"],ev(6),raw="99% ee",normalized=q(99,"%"),scope={"metric":"ee"}); outcome(r,"dr","selectivity",["m_2S_3S_3ha"],ev(6),raw="6:1 dr",normalized="6:1",scope={"metric":"dr"})
classify_partial(r,["C-C","N-H"],ev(6),"Partial net classification inferred from the coupling description; complete oxazoline/product mapping is unavailable.")

r=rxn("scheme6_deprotection_6a",(6,),categories=["derivatization"],name="protecting-group removal",transform="deprotection to alcohol 6a",location="page-6 prose introducing Scheme 6",
      unresolved=["Reagent, quantities, temperature, time, yield, workup, and purification are not stated on page 6."])
addp(r,"m_2S_3S_3ha","reactant",ev(6)); addp(r,"m_6a","product",ev(6),product_role="deprotected alcohol"); add_step(r,"reaction",ev(6),["p_m_2S_3S_3ha"],"remove the protecting group")
classify_partial(r,["O-H"],ev(6),"Only the prose-described conversion to an alcohol supports an O-H formation class; exact cleavage/mapping is unresolved.")

r=rxn("scheme6_metathesis_8a",(6,),categories=["derivatization"],name="olefin metathesis",transform="cross metathesis",location="page-6 prose introducing Scheme 6",
      unresolved=["Page 6 does not state charges, solvent, temperature, time, or yield measurement method."])
addp(r,"m_6a","reactant",ev(6)); addp(r,"m_7","reactant",ev(6),pid="p_7"); addp(r,"m_8a","product",ev(6),product_role="metathesis product"); addp(r,"m_HGII","catalyst",ev(6),pid="p_HGII")
add_step(r,"reaction",ev(6),["p_m_6a","p_7","p_HGII"],"olefin metathesis with Hoveyda-Grubbs II catalyst")
outcome(r,"yield","yield",["m_8a"],ev(6),raw="65%",normalized=q(65,"%")); classify_partial(r,["DB-C-C"],ev(6),"Net new alkene connectivity inferred from the stated olefin metathesis; complete mapping and coproducts are unavailable.")

r=rxn("scheme6_to_mycestericin_F",(6,),categories=["derivatization"],name="completion of mycestericin F synthesis",transform="hydrogenation and hydrolysis",location="page-6 prose introducing Scheme 6",
      unresolved=["The intermediate identity, reagents/conditions, stage yields, measurement methods, and exact hydrolysis sequence are not stated on page 6."])
addp(r,"m_8a","reactant",ev(6)); addp(r,"m_mycF","product",ev(6),product_role="natural product")
add_step(r,"hydrogenation",ev(6),["p_m_8a"],"hydrogenate the double bond"); add_step(r,"hydrolysis",ev(6),["p_m_8a"],"hydrolyze the oxazoline ring and ester group")
classify_partial(r,["C-H","O-H"],ev(6),"Only class-level hydrogenation/hydrolysis changes are inferred from the prose; exact atom mapping is unavailable.")

r=rxn("scheme6_syn_3hb",(6,),categories=["scale_up"],name="syn-selective coupling for mycestericin G",transform="alkoxyallene/oxazoline coupling",location="page-6 prose introducing Scheme 6",
      unresolved=["The source prints product label 3hb; relation to Table 3 numbering is unresolved. Yield, dr, ee, scale, and full conditions are not stated on page 6."])
addp(r,"m_1a","reactant",ev(6)); addp(r,"m_2h","reactant",ev(6),pid="p_2h"); addp(r,"m_3hb","product",ev(6),product_role="syn product")
addp(r,"m_Pd_L3","catalyst",ev(6)); addp(r,"m_L3","ligand",ev(6),pid="p_L3"); addp(r,"m_AgOAc","catalyst",ev(6),pid="p_Ag")
add_step(r,"reaction",ev(6),["p_m_1a","p_2h"],"couple using the L3-Pd/L3-Ag system")
classify_partial(r,["C-C","N-H"],ev(6),"Partial net classification inferred from the syn-coupling description; exact product identity/mapping is unresolved.")

r=rxn("scheme6_to_mycestericin_G",(6,),categories=["derivatization"],name="total synthesis of mycestericin G",transform="similar route from syn-coupling product",location="page-6 prose introducing Scheme 6",
      unresolved=["Page 6 states only that a similar route was used; intermediate identities, conditions, yields, workup, and exact product-label correspondence are unavailable."])
addp(r,"m_3hb","reactant",ev(6)); addp(r,"m_mycG","product",ev(6),product_role="natural product"); add_step(r,"route",ev(6),["p_m_3hb"],"convert by a route stated to be similar to the mycestericin F synthesis")
r["classification"]["coverage"]=F("insufficient transformation detail", "unresolved", ev(6), method="inference", derivation="The page supplies no unique atom-to-atom transformation for the similar route.")
r["classification"]["limitations"]=["No bond labels are asserted because the route is not chemically delimited on page 6."]


# Source context is deliberately separate from synthetic records.
bundle["source_context"] = [
    {"context_id":f"{DOC}:bibliographic_context","document_id":DOC,"kind":"bibliographic_context","evidence_ids":ev(1),
     "content":{"title":TITLE,"doi":"10.1002/anie.202014510","note":"Page 1 title/abstract/SI availability notice and illustrative mycestericin structures."}},
    {"context_id":f"{DOC}:structure_assignment_3ha","document_id":DOC,"kind":"structure_assignment","evidence_ids":ev(3,4),
     "content":{"compound":"(2S,3S)-3ha","method":"single-crystal X-ray diffraction","assignment":"(2S,3S)","ccdc":"2036502 is cited later; external CIF not supplied"}},
    {"context_id":f"{DOC}:mechanistic_discussion","document_id":DOC,"kind":"mechanistic_context","evidence_ids":ev(4,5,6),
     "content":{"note":"Author rationale involving Cu/Ag azomethine ylides and Pd-allyl conformers A-F; retained separately from executed experiments."}},
    {"context_id":f"{DOC}:negative_characterization_inventory_pages_1_6","document_id":DOC,"kind":"coverage_context","evidence_ids":ev(1,2,3,4,5,6),
     "content":{"note":"No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block occurs on PDF pages 1-6. Bibliographic and acknowledgment material is not a synthesis experiment."}},
]

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
(ROOT/"runs/current/submitted.json").write_text(json.dumps(bundle,ensure_ascii=False,indent=2)+"\n")

# Raw visual inventory and required coverage artifact.
page_regions = {
  1:[("title/abstract and introduction","bibliographic context; no concrete current-paper conditions"),("Scheme 1a","one grouped cited prior-work family"),("Scheme 1b and natural-product illustrations","current generic transformation; mycestericins are context")],
  2:[("Table 1 header and ligand structures","common substrates, catalysts, ligands, and conditions"),("Table 1 rows 1-22","22 independent optimization experiments"),("Table 1 footnotes","crude 1H NMR yields/dr, chiral-HPLC ee, isolated parenthetical yields, entry-20 loading change")],
  3:[("optimization discussion","entries 10-22 narrative and no-secondary-metal controls"),("Scheme 2","four stereoisomer outcomes; two merge with Table 1 entries 10 and 20"),("scope discussion","prose cross-check including the 3ja 56% versus Table 2 83% conflict and 3ha X-ray assignment")],
  4:[("Table 2 alkoxyallene scope","15 product arrows 3ab-3ap"),("Table 2 aldimine-ester scope","10 product arrows 3aa-3ja"),("Table 2 footnote","anti conditions, isolated yield, crude-NMR dr, chiral-HPLC ee"),("Scheme 4 discussion continuation","proposed stereochemical rationale only")],
  5:[("Table 3 alkoxyallene scope","13 product arrows 4ab-4an"),("Table 3 aldimine-ester scope","10 product arrows 4aa-4ja"),("Table 3 footnote","syn conditions, isolated yield, crude-NMR dr, chiral-HPLC ee"),("Scheme 3","eight stereoisomer outcomes; two merge with Table 3 4am/4an"),("mechanistic prose","continuation of proposed Scheme 4 rationale")],
  6:[("Scheme 4","proposed mechanism/rationalization, not an executed experiment"),("Scheme 5","five concrete derivatization arrows 5a-5e with conditions"),("mycestericin prose","six detail-limited concrete/route records preceding Scheme 6"),("conclusion/acknowledgments","non-experimental context")],
}

coverage={"document_id":DOC,"semantic_accuracy_verified":False,"pages":[]}
for page in range(1,7):
    rids=[r["reaction_id"] for r in reactions if I(page) in r["evidence_ids"]]
    unresolved=[]
    if page in (4,5): unresolved.append("Compact scope structures were not all converted to complete source-verified SMILES.")
    if page==6: unresolved.append("Scheme 6 is printed on PDF page 7 and is outside this session's claimed coverage; page-6 prose records retain missing conditions.")
    coverage["pages"].append({"page":page,"state":"attempted","regions":[{"region":a,"finding":b,"evidence_ids":[I(page),T(page)]} for a,b in page_regions[page]],"record_ids":rids,"unresolved":unresolved})
(ROOT/"runs/current/coverage/pages_0001_0006.json").write_text(json.dumps(coverage,ensure_ascii=False,indent=2)+"\n")

raw={"document_id":DOC,"page_range":[1,2,3,4,5,6],"visual_inspection":"original rendered page images opened at original detail",
     "inventories":{str(k):[{"region":a,"finding":b} for a,b in v] for k,v in page_regions.items()},
     "deduplication":["Scheme 2 standard products merged with Table 1 entries 10 and 20.","Scheme 3 standard L3-Pd/L3-Ag products 4am/4an merged with Table 3 rows.","Same label never merged across different conditions merely by product identity."],
     "conflicts":["Table 1 entry 10: 70% isolated in table versus 71% in prose/Scheme 2.","3ja: 83% isolated in Table 2 versus 56% in page-3 prose."],
     "outside_scope_context":"Other PDF pages were not claimed as covered; Scheme 6 on page 7 was not used to fill page-6 missing fields."}
(ROOT/"runs/current/raw_model/page_extraction.json").write_text(json.dumps(raw,ensure_ascii=False,indent=2)+"\n")

stereo=[]
for mid,m in molecules.items():
    smi=st.value(m["smiles"])
    if smi:
        stereo.append({"molecule_id":mid,"candidate_smiles":smi,"diagnostics":st.stereo_inspect(smi)})
(ROOT/"runs/current/raw_model/stereo_checks.json").write_text(json.dumps(stereo,ensure_ascii=False,indent=2)+"\n")

print(json.dumps({"molecules":len(bundle["molecules"]),"reactions":len(bundle["reactions"]),"source_context":len(bundle["source_context"])},ensure_ascii=False))
