#!/usr/bin/env python3
"""Build the evidence-linked v3.1 machine submission from pages 1--5.

This is an extraction artifact, not a chemistry validator.  Dense source-drawn
structures remain unresolved rather than being converted to guessed SMILES.
"""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_dcd68f4bd534e2f4"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 6)}
T = {p: f"{DOC}:page{p}:text" for p in range(1, 6)}


def q(number, unit, raw_text=None, qualifier=None):
    d = {"number": number, "unit": unit}
    if raw_text is not None:
        d["raw_text"] = raw_text
    if qualifier is not None:
        d["qualifier"] = qualifier
    return d


def fx(raw=None, normalized=None, ev=(), method="original_extraction", derivation=None,
       status=None):
    if normalized is None and raw is not None:
        normalized = raw
    return st.fact(raw, normalized, ev, method=method, derivation=derivation, status=status)


def inferred(value, ev, why):
    return fx(None, value, ev, method="inference", derivation=why)


def derived(value, ev, why):
    return fx(None, value, ev, method="structure_derivation", derivation=why)


bundle = st.read(RUN / "submitted.json")
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
molecules = {}
stereo_diagnostics = {}


def molecule(mid, label, ev, name=None, smiles=None, chirality="unresolved",
             geometry="unresolved", representation="unresolved", note=None,
             atom_centers=None, double_bonds=None, axial=None, mixture=None,
             assignment=None):
    m = st.new_molecule(mid, DOC, ev)
    m["occurrences"][0]["compound_label"] = fx(label, label, ev) if label is not None else fx()
    original = name or label
    m["names"]["original"] = fx(original, original, ev) if original else fx()
    m["names"]["english"] = fx(original, original, ev) if original else fx()
    m["names"]["abbreviations"] = fx()
    m["names"]["aliases"] = fx()
    m["cas_numbers"] = fx()
    m["identity_note"] = fx(note, note, ev) if note else fx()
    m["external_properties"] = []
    if smiles:
        diag = st.stereo_inspect(smiles)
        stereo_diagnostics[mid] = diag
        m["identity_status"] = "source_identified"
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smiles, ev,
            "Candidate SMILES was manually derived from the source name/drawing; RDKit parsing and stereochemistry diagnostics do not verify source identity.")
        m["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], ev, "Canonical isomeric SMILES calculated from the unverified candidate SMILES.")
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], ev, "Stereo graph key calculated from the unverified candidate SMILES.")
        if atom_centers is None:
            atom_centers = diag["atom_centers"] or None
        if double_bonds is None:
            double_bonds = diag["double_bonds"] or None
    else:
        m["identity_status"] = "source_identified"
        m["representation_type"] = "unresolved"
        m["smiles"] = fx()
        m["stereochemistry"]["isomeric_smiles"] = fx()
        m["stereochemistry"]["identity_key"] = fx()
    stereo_ev = ev
    def sf(value, why):
        return inferred(value, stereo_ev, why)
    m["stereochemistry"]["atom_centers"] = sf(atom_centers, "Source drawing/name and candidate graph were reviewed for atom-centered stereochemistry; any listed '?' remains unassigned.") if atom_centers is not None else fx()
    m["stereochemistry"]["double_bonds"] = sf(double_bonds, "Source drawing/name was reviewed for E/Z information; unspecified entries are not converted to E or Z.") if double_bonds is not None else fx()
    m["stereochemistry"]["axial"] = sf(axial, "Source drawing was reviewed for axial/allene stereochemistry.") if axial is not None else fx()
    m["stereochemistry"]["mixture"] = sf(mixture, "Source statement was reviewed for mixture/racemate information.") if mixture is not None else fx()
    m["stereochemistry"]["assignment"] = sf(assignment, "Source statement/drawing was reviewed for absolute or relative stereochemical assignment.") if assignment is not None else fx()
    m["stereochemistry"]["chirality_status"] = sf(chirality, "Stereochemistry status records source specification, not absence of stereo tags.")
    m["stereochemistry"]["geometry_status"] = sf(geometry, "Geometry status records source specification, not absence of slash notation.")
    m["stereochemistry"]["representation_status"] = sf(representation, "Representation status compares the source stereo information with the stored candidate representation.")
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


def unresolved(mid, label, ev, name=None, note=None, geometry="unspecified",
               double_bonds=None):
    return molecule(mid, label, ev, name=name, chirality="unspecified",
                    geometry=geometry, representation="not_encoded", note=note,
                    double_bonds=double_bonds)


# Core source-drawn molecules and common agents.
molecule("m_1", "1", [E[1], E[2]], "N-(o-iodophenyl)allenamide 1",
         "C=C=CN(c1ccccc1I)S(=O)(=O)c1ccc(C)cc1", chirality="achiral",
         geometry="not_applicable", representation="partial",
         axial=[{"axis": "allene C=C=C", "status": "not axially stereogenic", "basis": "one terminal allene carbon is CH2 in the source drawing"}],
         note="Candidate connectivity for source-drawn terminal allenamide 1; exact drawing remains machine-reviewable.")
molecule("m_1prime", "1′", [E[1], E[3]], "propargylamide 1′",
         "C#CCN(c1ccccc1I)S(=O)(=O)c1ccc(C)cc1", chirality="achiral",
         geometry="not_applicable", representation="encoded",
         note="Candidate connectivity for source-drawn propargylamide 1′.")
molecule("m_2", "2", [E[1], E[2]], "2-phenylpropanal",
         "CC(C=O)c1ccccc1", chirality="unspecified", geometry="not_applicable",
         representation="not_encoded", note="The source drawing has no wedge/dash or R/S assignment.")
unresolved("m_3", "3", [E[1], E[2], E[3]], "product 3",
           "Source-drawn N-tosyl indole aldehyde. CCDC 2046574 confirms a crystal structure, but the supplied pages do not print an absolute configuration; the achiral standard reaction and lack of stereo marks are not converted into a racemate claim.")

simple = [
    ("m_A1", "A1", "diphenylmethanamine", "NC(c1ccccc1)c1ccccc1", "achiral"),
    ("m_A2", "A2", "proline (configuration not specified)", "O=C(O)C1CCCN1", "unspecified"),
    ("m_A3", "A3", "pyrrolidine", "C1CCNC1", "achiral"),
    ("m_tBuONa", "t-BuONa", "sodium tert-butoxide", "[Na+].[O-]C(C)(C)C", "achiral"),
    ("m_Cs2CO3", "Cs2CO3", "cesium carbonate", "[Cs+].[Cs+].[O-]C([O-])=O", "achiral"),
    ("m_K3PO4", "K3PO4", "potassium phosphate", "[K+].[K+].[K+].[O-]P([O-])([O-])=O", "achiral"),
    ("m_Et3N", "Et3N", "triethylamine", "CCN(CC)CC", "achiral"),
    ("m_THF", "THF", "tetrahydrofuran", "C1CCOC1", "achiral"),
    ("m_PhMe", "PhMe", "toluene", "Cc1ccccc1", "achiral"),
    ("m_MeCN", "MeCN", "acetonitrile", "CC#N", "achiral"),
    ("m_DMSO", "DMSO", "dimethyl sulfoxide", "CS(C)=O", "achiral"),
    ("m_DBU", "DBU", "1,8-diazabicyclo[5.4.0]undec-7-ene", "C1CCC2=NCCCN2CC1", "achiral"),
    ("m_LDA", "LDA", "lithium diisopropylamide", "[Li+].[N-](C(C)C)C(C)C", "achiral"),
    ("m_MeOH", "MeOH", "methanol", "CO", "achiral"),
    ("m_PMBNH2", "PMBNH2", "(4-methoxyphenyl)methanamine", "COc1ccc(CN)cc1", "achiral"),
    ("m_Boc2O", "(Boc)2O", "di-tert-butyl dicarbonate", "CC(C)(C)OC(=O)OC(=O)OC(C)(C)C", "achiral"),
    ("m_K2CO3", "K2CO3", "potassium carbonate", "[K+].[K+].[O-]C([O-])=O", "achiral"),
    ("m_phosphonate", "(EtO)2POCH2COOEt", "triethyl phosphonoacetate", "CCOC(=O)CP(=O)(OCC)OCC", "achiral"),
]
for mid, label, name, smi, ch in simple:
    ev = [E[1]] if mid in {"m_A1","m_A2","m_A3","m_tBuONa","m_Cs2CO3","m_K3PO4","m_Et3N","m_THF","m_PhMe","m_MeCN","m_DMSO"} else [E[2], E[3]]
    molecule(mid, label, ev, name, smi, chirality=ch, geometry="not_applicable",
             representation="encoded" if ch == "achiral" else "not_encoded")

for mid, label, name, ev in [
    ("m_PdPPh3", "Pd(PPh3)4", "tetrakis(triphenylphosphine)palladium(0)", [E[1],E[2]]),
    ("m_PdPtBu3", "Pd(PtBu3)2", "bis(tri-tert-butylphosphine)palladium(0)", [E[1]]),
    ("m_Pd2dba3", "Pd2(dba)3", "tris(dibenzylideneacetone)dipalladium", [E[1]]),
    ("m_Ar", "Ar", "argon", [E[1],E[2]]),
    ("m_NaBH3CN", "NaBH3CN", "sodium cyanoborohydride", [E[3]]),
    ("m_MeMgBr", "MeMgBr", "methylmagnesium bromide", [E[3]]),
    ("m_NiCl2_6H2O", "NiCl2·6H2O", "nickel(II) chloride hexahydrate", [E[3]]),
    ("m_NaBH4", "NaBH4", "sodium borohydride", [E[3]]),
    ("m_unknown_chiral_ligand", "chiral ligand", "unidentified chiral ligand from screening", [E[3]]),
    ("m_unknown_chiral_amine", "chiral amine catalyst", "unidentified chiral amine catalyst from screening", [E[3]]),
]:
    molecule(mid, label, ev, name=name, note="Identity/connectivity is not fully specified on the supplied pages; no SI was associated.")

# Source-drawn scope inputs and products.  Distinct drawings stay distinct even
# when the source does not print input compound numbers.
for n in range(4, 48):
    unresolved(f"m_{n}", str(n), [E[2], E[3]] if n in (39,) else [E[2]], f"product {n}",
               "Complete source drawing is retained in page-image evidence. No absolute/relative configuration is printed; no complete SMILES was guessed.")
for n in range(48, 52):
    unresolved(f"m_{n}", str(n), [E[3]], f"product {n}",
               "Complete source drawing is retained in Scheme 5. No absolute/relative configuration is printed; no complete SMILES was guessed.",
               geometry="unspecified" if n == 49 else "not_applicable",
               double_bonds=[{"feature": "new alkene in product 49", "configuration": "unspecified by source"}] if n == 49 else None)

for n in list(range(4, 11)) + list(range(26, 36)):
    unresolved(f"m_allenamide_for_{n}", None, [E[2]], f"source-drawn allenamide leading to product {n}",
               "Unnumbered concrete allenamide drawing in the scope; retained as a distinct unresolved molecule.")
for n in list(range(11, 26)) + list(range(27, 30)):
    unresolved(f"m_aldehyde_for_{n}", None, [E[2]], f"source-drawn aldehyde leading to product {n}",
               "Unnumbered concrete aldehyde drawing inferred from the corresponding product row; no SMILES guessed.")
for n in range(36, 48):
    unresolved(f"m_pronucleophile_for_{n}", None, [E[2]], f"source-drawn pronucleophile leading to product {n}",
               "Unnumbered concrete pronucleophile drawing inferred from the corresponding product row; no SMILES guessed.")


def participant(pid, mid, role, ev, amount=None, equivalents=None, mol_percent=None,
                volume=None, ref=None, product_role=None, stage="reaction"):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    p["product_role"] = fx(product_role, product_role, ev) if product_role else fx()
    if amount is not None:
        p["charges"]["amount"] = fx(f"{amount} mmol", q(amount,"mmol"), ev)
    if equivalents is not None:
        p["charges"]["equivalents"] = fx(f"{equivalents} equiv", q(equivalents,"equiv"), ev)
    if mol_percent is not None:
        p["charges"]["mol_percent"] = fx(f"{mol_percent} mol %", q(mol_percent,"mol%"), ev)
    if volume is not None:
        p["charges"]["volume"] = fx(f"{volume} mL", q(volume,"mL"), ev)
    if ref is not None:
        p["equivalence_reference"] = inferred(ref, ev,
            "The paper states the loading/equivalents but does not print the denominator; the reference follows the reaction-scale substrate convention and remains machine-reviewable.")
    return p


def change(change_id, a, b, before, after, kind):
    return {"change_id": change_id, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
            "before_order": before, "after_order": after, "change_type": kind}


def set_classification(r, kind, ev):
    c = r["classification"]
    c["mapped_reaction_smiles"] = fx()
    limitations = ["No atom-complete mapping was asserted for dense source drawings.",
                   "Labels are category-level structure derivations from the depicted net transformation and are partial."]
    if kind == "none":
        c["bond_changes"] = []
        c["labels"] = inferred([], ev, "No observed product transformation is classified for this no-reaction/proposed-mechanism record.")
        c["coverage"] = inferred("not_applicable", ev, "No observed net bond formation is applicable.")
    else:
        specs = {
            "cascade": [change("b1","C","C",0,1,"formed"), change("b2","C","C",0,1,"formed")],
            "benzylation": [change("b1","C","C",0,1,"formed")],
            "reductive_amination": [change("b1","C","N",0,1,"formed"), change("b2","C","H",0,1,"hydrogen_gain")],
            "olefination": [change("b1","C","C",0,2,"formed")],
            "detosylation": [change("b1","N","H",0,1,"hydrogen_gain")],
            "nitrile_reduction_boc": [change("b1","C","N",0,1,"formed"), change("b2","C","H",0,1,"hydrogen_gain"), change("b3","C","N",3,1,"order_decreased")],
        }[kind]
        c["bond_changes"] = [derived(x, ev, "Category-level bond change inferred from the source-drawn reactant/product structures; atom identities are unresolved.") for x in specs]
        labels = {}
        for x in specs:
            if x["change_type"] not in ("formed","hydrogen_gain","order_increased"):
                continue
            lab = st.bond_label(x["atoms"][0]["element"], x["atoms"][1]["element"], x["after_order"])
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(x["change_id"])
        c["labels"] = derived(list(labels.values()), ev, "Taxonomy labels calculated from the explicitly stored partial bond changes.")
        c["coverage"] = inferred("partial", ev, "Only unambiguous category-level changes were encoded; hydrogen/atom provenance and full mapping remain unresolved.")
    c["limitations"] = limitations
    c["model_votes"] = []


def reaction(rid, record_type, categories, ev, location, name, transformation, class_kind,
             unresolved_items):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = fx(categories, categories, ev)
    r["location"] = fx(location, location, ev)
    r["reaction_name"] = fx(name, name, ev) if name else fx()
    r["transformation_type"] = inferred(transformation, ev, "Transformation type normalized from the source title, scheme caption, and depicted structures.") if transformation else fx()
    r["general_scheme_ids"] = []
    r["reaction_scale"] = fx()
    r["solvent_conditions"] = {k: fx() for k in ("composition","mixture_ratio","total_volume","concentration","concentration_basis")}
    r["special_conditions"] = {k: fx() for k in ("light_source","wavelength","power","current","voltage","electrodes","flow_parameters","reactor")}
    r["procedure"] = {k: fx() for k in ("original_text","general_procedure_reference","modifications","addition_order","pretreatment")}
    r["participants"] = []
    r["steps"] = []
    r["outcomes"] = []
    r["notes"] = []
    r["review"] = {"status":"needs_review","method":"model","reviewer":None,"checked_fields":[],"unresolved":list(unresolved_items)}
    r["completeness"] = "partial"
    set_classification(r, class_kind, ev)
    bundle["reactions"].append(r)
    return r


def standard_step(r, ev, pids, temp=50, time=8, atmosphere="Ar", operation=None):
    s = st.new_step("reaction", ev)
    s["participant_ids"] = pids
    s["operation"] = fx(operation, operation, ev) if operation else fx()
    s["temperature"] = fx(f"{temp} °C", q(temp,"°C"), ev) if temp is not None else fx()
    s["time"] = fx(f"{time} h", q(time,"h"), ev) if time is not None else fx()
    s["atmosphere"] = fx(atmosphere, atmosphere, ev) if atmosphere else fx()
    s["pressure"] = fx()
    r["steps"].append(s)


def outcome(r, oid, kind, targets, raw, normalized, ev, basis=None, yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["value"] = fx(raw, normalized, ev)
    o["product_role"] = fx("main product", "main_product", ev) if targets else fx()
    o["yield_type"] = fx(yield_type, yield_type, ev) if yield_type else fx()
    o["measurement_basis"] = fx(basis, basis, ev) if basis else fx()
    o["scope"] = fx(scope, scope, ev) if scope else fx()
    r["outcomes"].append(o)


# Scheme 1: one cited generic example, two current generic transformations, and
# an author-proposed pathway.  None are counted as concrete present-paper runs.
r = reaction("scheme1_previous_benzylation", "literature_example", ["other"], [E[1]],
             {"page":1,"region":"Scheme 1a"}, "benzylation of active methine compounds",
             "generic benzylation", "benzylation", ["Generic R groups and several leaving groups/catalyst classes prevent a concrete experiment identity."])
r["notes"].append(fx("Previous works", "Scheme explicitly marks this as previous work, not a present-paper experiment.", [E[1]]))

for suffix, title in [("indole","allenamide carbopalladation/allylation to indoles"),
                      ("isoquinolinone","allenamide carbopalladation/allylation to isoquinolinones")]:
    r = reaction(f"scheme1_current_{suffix}", "generic_scheme", ["other"], [E[1],T[2]],
                 {"page":1,"region":"Scheme 1b"}, title,
                 "cascade oxidative addition/carbopalladation/allylation", "cascade",
                 ["Generic Ar/R/EWG variables and PdL intermediates are not concrete molecule identities."])

r = reaction("proposed_pathway", "proposed_mechanism", ["other"], [E[1],T[2],T[3]],
             {"page":1,"region":"Scheme 1b intermediates; prose pages 2-3"},
             "sequential oxidative addition, carbopalladation, and allylic alkylation",
             "author-proposed cascade pathway", "none",
             ["This is an author-proposed sequence, not a separately demonstrated mechanistic experiment."])


# Table 1, all 14 rows.  Product 3 in Scheme 2 at 92% is the same standard run
# as entry 1 and contributes page-2 evidence rather than a duplicate record.
opt_rows = {
  1:("none",92,"yield"), 2:("Pd(PtBu3)2 instead of Pd(PPh3)4",38,"yield"),
  3:("Pd2(dba)3 instead of Pd(PPh3)4","trace","trace"), 4:("A2 instead of A1",21,"yield"),
  5:("A3 instead of A1",26,"yield"), 6:("Cs2CO3 instead of t-BuONa",51,"yield"),
  7:("K3PO4 instead of t-BuONa",12,"yield"), 8:("Et3N instead of t-BuONa","trace","trace"),
  9:("PhMe instead of THF",36,"yield"), 10:("MeCN instead of THF",45,"yield"),
 11:("DMSO instead of THF",23,"yield"), 12:("without Pd(PPh3)4","NR","no_reaction"),
 13:("without A1",40,"yield"), 14:("1′ instead of 1",53,"yield")}

for i,(deviation,val,kind) in opt_rows.items():
    ev = [E[1],T[2],E[2]] if i == 1 else [E[1],T[2]] if i in range(2,14) else [E[1],T[3]]
    r = reaction(f"table1_entry_{i:02d}", "specific_experiment", ["optimization"] + (["control"] if i in (12,13) else []), ev,
                 {"page":1,"region":f"Table 1 entry {i}"}, "Optimization of Reaction Conditions",
                 "allenamide carbopalladation/allylation", "none" if i == 12 else "cascade",
                 ["SI experimental procedure/workup/purification is unavailable."])
    r["general_scheme_ids"] = [f"{DOC}:scheme1_current_indole"]
    substrate = "m_1prime" if i == 14 else "m_1"
    r["participants"] += [participant("p_allenamide",substrate,"reactant",ev,amount=0.2), participant("p_aldehyde","m_2","reactant",ev,amount=0.4)]
    catalyst = {2:"m_PdPtBu3",3:"m_Pd2dba3"}.get(i,"m_PdPPh3")
    if i != 12:
        r["participants"].append(participant("p_palladium",catalyst,"catalyst",ev,mol_percent=10,ref=substrate))
    amine = {4:"m_A2",5:"m_A3"}.get(i,"m_A1")
    if i != 13:
        r["participants"].append(participant("p_amine",amine,"catalyst",ev,mol_percent=30,ref=substrate))
    base = {6:"m_Cs2CO3",7:"m_K3PO4",8:"m_Et3N"}.get(i,"m_tBuONa")
    r["participants"].append(participant("p_base",base,"reagent",ev,amount=0.2))
    solvent = {9:"m_PhMe",10:"m_MeCN",11:"m_DMSO"}.get(i,"m_THF")
    r["participants"].append(participant("p_solvent",solvent,"solvent",ev,volume=1.0))
    if i != 12:
        r["participants"].append(participant("p_product","m_3","product",ev,product_role="desired product"))
    r["reaction_scale"] = fx("1 (0.2 mmol)", {"reference_molecule_id":substrate,"amount":q(0.2,"mmol")}, ev)
    r["solvent_conditions"]["composition"] = fx(molecules[solvent]["names"]["original"]["normalized_value"], solvent, ev)
    r["solvent_conditions"]["total_volume"] = fx("1.0 mL", q(1.0,"mL"), ev)
    r["special_conditions"]["reactor"] = fx("sealed tube", "sealed tube", ev)
    r["procedure"]["general_procedure_reference"] = fx("standard conditions; deviation: "+deviation, deviation, ev)
    standard_step(r,ev,[p["participant_id"] for p in r["participants"] if p["stage"]=="reaction"],operation="sealed-tube reaction")
    if kind == "yield":
        outcome(r,"yield","yield",["m_3"],f"{val}%",q(val,"%"),ev,basis="isolated",yield_type="isolated")
    elif kind == "trace":
        outcome(r,"trace","trace",["m_3"],"trace","trace",ev,basis="Table 1 footnote: isolated yields")
    else:
        outcome(r,"no_reaction","no_reaction",[],"NR; no reaction","no reaction",ev,basis="NR, no reaction")
    r["notes"].append(fx(deviation, {"deviation_from_standard":deviation}, ev))


scope2_yields = {4:98,5:91,6:92,7:78,8:94,9:96,10:85,11:94,12:98,13:95,14:75,15:95,16:90,17:67,18:92,19:87,20:94,21:97,22:96,23:91,24:98,25:98}
for n,y in scope2_yields.items():
    r = reaction(f"scheme2_product_{n}","specific_experiment",["substrate_scope"],[E[2],T[3]],
                 {"page":2,"region":f"Scheme 2 product {n}"},"Substrate Scope for the Synthesis of Indoles",
                 "allenamide carbopalladation/allylation","cascade",["SI workup/purification and compound characterization are unavailable.","Dense structures lack atom-complete machine encoding."])
    r["general_scheme_ids"]=[f"{DOC}:scheme1_current_indole"]
    a = f"m_allenamide_for_{n}" if n <= 10 else "m_1"
    ald = "m_2" if n <= 10 else f"m_aldehyde_for_{n}"
    r["participants"]=[participant("p_allenamide",a,"reactant",[E[2]],amount=0.2),participant("p_aldehyde",ald,"reactant",[E[2]],amount=0.4),participant("p_pd","m_PdPPh3","catalyst",[E[2]],mol_percent=10,ref=a),participant("p_A1","m_A1","catalyst",[E[2]],mol_percent=30,ref=a),participant("p_base","m_tBuONa","reagent",[E[2]],amount=0.2),participant("p_solvent","m_THF","solvent",[E[2]],volume=1.0),participant("p_product",f"m_{n}","product",[E[2]],product_role="main product")]
    r["reaction_scale"]=fx("allenamide (0.2 mmol)",{"reference_molecule_id":a,"amount":q(0.2,"mmol")},[E[2]])
    r["solvent_conditions"]["composition"]=fx("THF","m_THF",[E[2]])
    r["solvent_conditions"]["total_volume"]=fx("1.0 mL",q(1.0,"mL"),[E[2]])
    r["procedure"]["general_procedure_reference"]=fx("Scheme 2 footnote a", "Scheme 2 general conditions", [E[2]])
    standard_step(r,[E[2]],[p["participant_id"] for p in r["participants"]])
    outcome(r,"yield","yield",[f"m_{n}"],f"{y}%",q(y,"%"),[E[2]],basis="isolated",yield_type="isolated")

# Scheme 3: product-specific paired inputs retained unresolved.
scope3_yields={26:96,27:87,28:87,29:84,30:91,31:98,32:57,33:58,34:64,35:90}
for n,y in scope3_yields.items():
    a=f"m_allenamide_for_{n}"
    ald="m_2" if n in (26,30,31,32,33,34,35) else f"m_aldehyde_for_{n}"
    r=reaction(f"scheme3_product_{n}","specific_experiment",["substrate_scope"],[E[2],T[3]],{"page":2,"region":f"Scheme 3 product {n}"},"Substrate Scope for the Synthesis of Isoquinolinones","allenamide carbopalladation/allylation","cascade",["SI workup/purification and compound characterization are unavailable.","Dense structures lack atom-complete machine encoding."])
    r["general_scheme_ids"]=[f"{DOC}:scheme1_current_isoquinolinone"]
    r["participants"]=[participant("p_allenamide",a,"reactant",[E[2]],amount=0.2),participant("p_aldehyde",ald,"reactant",[E[2]],amount=0.4),participant("p_pd","m_PdPPh3","catalyst",[E[2]],mol_percent=10,ref=a),participant("p_A1","m_A1","catalyst",[E[2]],mol_percent=30,ref=a),participant("p_base","m_tBuONa","reagent",[E[2]],amount=0.2),participant("p_solvent","m_THF","solvent",[E[2]],volume=1.0),participant("p_product",f"m_{n}","product",[E[2]],product_role="main product")]
    r["reaction_scale"]=fx("allenamide (0.2 mmol)",{"reference_molecule_id":a,"amount":q(0.2,"mmol")},[E[2]])
    r["solvent_conditions"]["composition"]=fx("THF","m_THF",[E[2]])
    r["solvent_conditions"]["total_volume"]=fx("1.0 mL",q(1.0,"mL"),[E[2]])
    r["procedure"]["general_procedure_reference"]=fx("Scheme 3 footnote a","Scheme 3 general conditions",[E[2]])
    standard_step(r,[E[2]],[p["participant_id"] for p in r["participants"]])
    outcome(r,"yield","yield",[f"m_{n}"],f"{y}%",q(y,"%"),[E[2]],basis="isolated",yield_type="isolated")

# Scheme 4, preserving every lettered exception in structured participants.
scope4_yields={36:85,37:54,38:77,39:79,40:82,41:43,42:44,43:93,44:55,45:68,46:67,47:60}
for n,y in scope4_yields.items():
    nuc=f"m_pronucleophile_for_{n}"
    r=reaction(f"scheme4_product_{n}","specific_experiment",["substrate_scope"],[E[2],T[3]],{"page":2,"region":f"Scheme 4 product {n} and footnotes"},"Substrate Scope of Pronucleophiles","allenamide carbopalladation/allylation","cascade",["SI workup/purification and compound characterization are unavailable.","Dense structures lack atom-complete machine encoding."])
    r["general_scheme_ids"]=[f"{DOC}:scheme1_current_indole"]
    allen_eq=2.0 if n in (46,47) else None
    nuc_eq=4.0 if n==45 else None
    r["participants"]=[participant("p_allenamide","m_1","reactant",[E[2]],amount=None if allen_eq else 0.2,equivalents=allen_eq,ref=nuc if allen_eq else None),participant("p_pronucleophile",nuc,"reactant",[E[2]],amount=None if nuc_eq else 0.4,equivalents=nuc_eq,ref="m_1" if nuc_eq else None),participant("p_pd","m_PdPPh3","catalyst",[E[2]],mol_percent=10,ref="m_1")]
    base="m_DBU" if n==38 else "m_LDA" if n in (41,42) else "m_Cs2CO3"
    if n==45:
        r["participants"].append(participant("p_base",base,"reagent",[E[2]],equivalents=1.0,ref="m_1"))
    else:
        r["participants"].append(participant("p_base",base,"reagent",[E[2]],amount=0.4))
    if n==47:
        r["participants"].append(participant("p_A1","m_A1","catalyst",[E[2]],mol_percent=30,ref="m_1"))
        r["review"]["unresolved"].append("The denominator for A1 30 mol% is not explicit in Scheme 4 footnote h.")
    r["participants"] += [participant("p_solvent","m_THF","solvent",[E[2]],volume=1.0),participant("p_product",f"m_{n}","product",[E[2]],product_role="main product")]
    scale_ref=nuc if n in (46,47) else "m_1"
    scale_raw="pronucleophile (0.4 mmol); allenamide 1 (2.0 equiv)" if n in (46,47) else "allenamide 1 (0.2 mmol)"
    r["reaction_scale"]=fx(scale_raw,{"reference_molecule_id":scale_ref,"amount":q(0.4 if n in (46,47) else 0.2,"mmol")},[E[2]])
    r["solvent_conditions"]["composition"]=fx("THF","m_THF",[E[2]])
    r["solvent_conditions"]["total_volume"]=fx("1.0 mL",q(1.0,"mL"),[E[2]])
    exceptions=[]
    if n==38: exceptions.append("c DBU (0.4 mmol)")
    if n in (41,42): exceptions.append("d LDA (0.4 mmol)")
    if n==45: exceptions += ["e Cs2CO3 (1.0 equiv)","f pronucleophile (4.0 equiv)"]
    if n in (46,47): exceptions.append("g allenamide 1 (2.0 equiv)")
    if n==47: exceptions.append("h A1 (30 mol %)")
    r["procedure"]["general_procedure_reference"]=fx("Scheme 4 footnote a"+("; "+"; ".join(exceptions) if exceptions else ""),{"general":"Scheme 4 conditions","exceptions":exceptions},[E[2]])
    standard_step(r,[E[2]],[p["participant_id"] for p in r["participants"]])
    outcome(r,"yield","yield",[f"m_{n}"],f"{y}%",q(y,"%"),[E[2]],basis="isolated",yield_type="isolated")


# Two distinct scale-up experiments embedded in scope schemes.
for n,scale,mass,y in [(3,2.5,1.02,98),(39,1.0,0.34,83)]:
    ev=[E[2],T[3]]
    r=reaction(f"scaleup_product_{n}","specific_experiment",["scale_up"],ev,{"page":2,"region":f"Scheme {'2' if n==3 else '4'} scale-up annotation for product {n}"},"Scale-up under optimal conditions","scale-up allenamide carbopalladation/allylation","cascade",["Only scale, product mass, yield, and reference to optimal conditions are supplied; SI procedure/workup is unavailable."])
    r["general_scheme_ids"]=[f"{DOC}:scheme1_current_indole"]
    second="m_2" if n==3 else "m_pronucleophile_for_39"
    r["participants"]=[participant("p_allenamide","m_1","reactant",ev,amount=scale),participant("p_partner",second,"reactant",ev),participant("p_product",f"m_{n}","product",ev,product_role="main product")]
    r["reaction_scale"]=fx(f"{scale} mmol",{"reference_molecule_id":"m_1","amount":q(scale,"mmol")},ev)
    r["procedure"]["general_procedure_reference"]=fx("under the optimal conditions","optimal conditions",ev)
    standard_step(r,ev,["p_allenamide","p_partner","p_product"],operation="optimal conditions; exact scale-up charges not fully printed")
    outcome(r,"product_amount","product_amount",[f"m_{n}"],f"{mass} g",q(mass,"g"),ev,basis="isolated material")
    outcome(r,"yield","yield",[f"m_{n}"],f"{y}%",q(y,"%"),ev,basis="isolated",yield_type="isolated")


# Scheme 5: four arrows, each a separate experiment.
r=reaction("scheme5_product_48","specific_experiment",["derivatization"],[E[3],T[3]],{"page":3,"region":"Scheme 5 first arrow"},"condensation followed by reduction","one-pot reductive amination","reductive_amination",["Reagent charges, stage durations, workup, and purification are unavailable."])
r["participants"]=[participant("p_3","m_3","reactant",[E[3]]),participant("p_PMBNH2","m_PMBNH2","reagent",[E[3]]),participant("p_MeOH","m_MeOH","solvent",[E[3]]),participant("p_NaBH3CN","m_NaBH3CN","reagent",[E[3]]),participant("p_48","m_48","product",[E[3]],product_role="main product")]
s1=st.new_step("condensation",[E[3]]);s1["participant_ids"]=["p_3","p_PMBNH2","p_MeOH"];s1["operation"]=fx("PMBNH2, MeOH, 100 °C","condensation",[E[3]]);s1["temperature"]=fx("100 °C",q(100,"°C"),[E[3]]);s1["time"]=fx();s1["atmosphere"]=fx();s1["pressure"]=fx();r["steps"].append(s1)
s2=st.new_step("reduction",[E[3]]);s2["participant_ids"]=["p_NaBH3CN","p_48"];s2["operation"]=fx("then NaBH3CN, 0 °C","NaBH3CN reduction",[E[3]]);s2["temperature"]=fx("0 °C",q(0,"°C"),[E[3]]);s2["time"]=fx();s2["atmosphere"]=fx();s2["pressure"]=fx();r["steps"].append(s2)
outcome(r,"yield","yield",["m_48"],"54%",q(54,"%"),[E[3],T[3]],basis="isolated",yield_type="isolated one-pot yield",scope="from product 3 to product 48 in one pot")

r=reaction("scheme5_product_49","specific_experiment",["derivatization"],[E[3],T[3]],{"page":3,"region":"Scheme 5 second arrow"},"Wittig reaction","olefination","olefination",["Reagent charges, reflux temperature/duration, E/Z geometry, workup, purification, and yield measurement basis are unavailable."])
r["participants"]=[participant("p_3","m_3","reactant",[E[3]]),participant("p_MeMgBr","m_MeMgBr","reagent",[E[3]]),participant("p_phosphonate","m_phosphonate","reagent",[E[3]]),participant("p_THF","m_THF","solvent",[E[3]]),participant("p_49","m_49","product",[E[3]],product_role="main product")]
standard_step(r,[E[3]],[p["participant_id"] for p in r["participants"]],temp=None,time=None,atmosphere=None,operation="MeMgBr and (EtO)2POCH2COOEt in THF at reflux")
r["steps"][0]["temperature"]=fx("reflux","reflux",[E[3]])
outcome(r,"yield","yield",["m_49"],"88%",q(88,"%"),[E[3],T[3]],scope="product 3 to product 49")

r=reaction("scheme5_product_50","specific_experiment",["derivatization"],[E[3],T[3]],{"page":3,"region":"Scheme 5 third arrow"},"N-tosyl deprotection","detosylation","detosylation",["Reagent charges, duration, workup, purification, and yield measurement basis are unavailable."])
r["participants"]=[participant("p_3","m_3","reactant",[E[3]]),participant("p_K2CO3","m_K2CO3","reagent",[E[3]]),participant("p_MeOH","m_MeOH","solvent",[E[3]]),participant("p_50","m_50","product",[E[3]],product_role="main product")]
standard_step(r,[E[3]],[p["participant_id"] for p in r["participants"]],temp=80,time=None,atmosphere=None,operation="K2CO3 in MeOH")
outcome(r,"yield","yield",["m_50"],"62%",q(62,"%"),[E[3],T[3]],scope="product 3 to product 50")

r=reaction("scheme5_product_51","specific_experiment",["derivatization"],[E[3],T[3]],{"page":3,"region":"Scheme 5 fourth arrow"},"nitrile reduction and Boc protection","reduction/protection to beta-amino acid derivative","nitrile_reduction_boc",["Reagent charges, duration, workup, purification, and yield measurement basis are unavailable."])
r["participants"]=[participant("p_39","m_39","reactant",[E[3]]),participant("p_Ni","m_NiCl2_6H2O","reagent",[E[3]]),participant("p_Boc","m_Boc2O","reagent",[E[3]]),participant("p_NaBH4","m_NaBH4","reagent",[E[3]]),participant("p_MeOH","m_MeOH","solvent",[E[3]]),participant("p_51","m_51","product",[E[3]],product_role="main product")]
standard_step(r,[E[3]],[p["participant_id"] for p in r["participants"]],temp=0,time=None,atmosphere=None,operation="NiCl2·6H2O, (Boc)2O, NaBH4 in MeOH; initial low-temperature stage")
s2=st.new_step("warm_to_room_temperature",[E[3]]);s2["participant_ids"]=[p["participant_id"] for p in r["participants"]];s2["operation"]=fx("0 °C - rt","warm from 0 °C to room temperature; transition time unreported",[E[3]]);s2["temperature"]=fx("rt","room temperature",[E[3]]);s2["time"]=fx();s2["atmosphere"]=fx();s2["pressure"]=fx();r["steps"].append(s2)
outcome(r,"yield","yield",["m_51"],"63%",q(63,"%"),[E[3],T[3]],scope="product 39 to beta-amino acid derivative 51")


# The main text reports one best asymmetric result but sends all identities and
# screening details to unavailable SI.  Preserve the concrete outcomes without
# inventing a ligand/catalyst or absolute configuration.
r=reaction("asymmetric_summary_product_3","specific_experiment",["other"],[T[3],E[3]],{"page":3,"region":"paragraph beside Scheme 5"},"asymmetric version of this method","asymmetric allenamide carbopalladation/allylation","cascade",["Exact chiral ligand, chiral amine catalyst, quantities, temperature, time, workup, and purification are referred to unavailable SI.","The major enantiomer and absolute configuration are not assigned on the supplied pages.","Yield measurement basis and ee measurement method are unreported."])
r["general_scheme_ids"]=[f"{DOC}:scheme1_current_indole"]
r["participants"]=[participant("p_1","m_1","reactant",[T[3]]),participant("p_2","m_2","reactant",[T[3]]),participant("p_ligand","m_unknown_chiral_ligand","ligand",[T[3]]),participant("p_amine","m_unknown_chiral_amine","catalyst",[T[3]]),participant("p_3","m_3","product",[T[3]],product_role="enantioenriched product; absolute configuration unassigned")]
standard_step(r,[T[3],E[3]],[p["participant_id"] for p in r["participants"]],temp=None,time=None,atmosphere=None,operation="asymmetric version after screening chiral ligands and chiral amine catalysts; exact best conditions in unavailable SI")
outcome(r,"yield","yield",["m_3"],"89%",q(89,"%"),[T[3]],scope="asymmetric version product 3")
outcome(r,"ee","selectivity",["m_3"],"40% ee",q(40,"%"),[T[3]],scope={"metric":"ee","major_enantiomer":"unassigned","order":"unreported"})


# Non-synthetic/contextual observations and explicit source limitations.
bundle["source_context"] = [
  {"context_id":f"{DOC}:graphical_abstract","kind":"graphical_summary","evidence_ids":[E[1]],"note":"Duplicates the generic present-work transformations; not a separate experiment."},
  {"context_id":f"{DOC}:xray_3","kind":"characterization","evidence_ids":[E[2],T[3]],"note":"Product 3 structure confirmed by single-crystal X-ray diffraction; CCDC 2046574. No standalone synthesis record created."},
  {"context_id":f"{DOC}:si_notice","kind":"source_limitation","evidence_ids":[E[3],T[3]],"note":"The article states SI exists, but no SI file was supplied or explicitly associated in this workspace."},
  {"context_id":f"{DOC}:references_pages_4_5","kind":"bibliography","evidence_ids":[E[4],T[4],E[5],T[5]],"note":"Acknowledgments and references only; cited reactions are not imported as present-paper experiments."},
]

# Raw transcription remains separate from the normalized bundle.
raw = {
  "document_id":DOC,
  "pages_reviewed":[1,2,3,4,5],
  "table1": [{"entry":i,"deviation":d,"reported_result":v,"result_kind":k} for i,(d,v,k) in opt_rows.items()],
  "scheme2_yields":scope2_yields,
  "scheme3_yields":scope3_yields,
  "scheme4_yields":scope4_yields,
  "scheme4_footnotes":{"c":"DBU (0.4 mmol)","d":"LDA (0.4 mmol)","e":"Cs2CO3 (1.0 equiv)","f":"Pronucleophile (4.0 equiv)","g":"Allenamide 1 (2.0 equiv)","h":"A1 (30 mol %)"},
  "scaleups":{"3":{"mass":"1.02 g","yield":"98%","scale":"2.5 mmol"},"39":{"mass":"0.34 g","yield":"83%","scale":"1.0 mmol"}},
  "scheme5":{"48":{"yield":"54%","basis":"isolated in one pot"},"49":{"yield":"88%","basis":None},"50":{"yield":"62%","basis":None},"51":{"yield":"63%","basis":None}},
  "asymmetric_summary":{"product":"3","yield":"89%","ee":"40% ee","major_enantiomer":None,"method":None,"conditions":"see unavailable SI"},
  "limitations":["No SI supplied/associated","Dense structures retained without guessed complete SMILES","Machine records remain needs_review"]
}
st.write(RUN / "raw_model/transcription.json", raw)
st.write(RUN / "raw_model/stereo_diagnostics.json", stereo_diagnostics)
st.write(RUN / "raw_model/source_mentions.json", st.source_mentions(
    {"documents": bundle["documents"]},
    ["1′", "product 3", "compound 39", "product 39", "40% ee", "CCDC 2046574"],
    context=300))

coverage_regions = {
  1: [
    {"region":"title/abstract/graphical abstract","finding":"Metadata and graphical summary; the graphical abstract duplicates the two present-work generic transformations and is not a separate experiment.","evidence_ids":[E[1],T[1]]},
    {"region":"Scheme 1a","finding":"One generic previous-work benzylation example; classified as literature_example, not as a present-paper experiment.","evidence_ids":[E[1],T[1]]},
    {"region":"Scheme 1b","finding":"Two current generic cascades (indole and isoquinolinone) plus schematic Pd intermediates supporting an author-proposed pathway.","evidence_ids":[E[1],T[1]]},
    {"region":"Table 1 entries 1-14 and footnotes","finding":"Fourteen independent optimization/control rows transcribed; entry 1 is the same 92% standard experiment repeated as product 3 in Scheme 2.","evidence_ids":[E[1],T[1]]},
  ],
  2: [
    {"region":"Scheme 2 products 3-25 and footnotes","finding":"Standard product 3 duplicate, distinct 2.5 mmol scale-up, and 22 independent product 4-25 substrate-scope experiments; all reported yields are isolated.","evidence_ids":[E[2],T[2]]},
    {"region":"Scheme 3 products 26-35 and footnotes","finding":"Ten independent isoquinolinone substrate-scope experiments with isolated yields.","evidence_ids":[E[2],T[2]]},
    {"region":"Scheme 4 products 36-47 and footnotes c-h","finding":"Twelve independent pronucleophile-scope experiments, all lettered condition exceptions structured; product 39 also has a distinct 1.0 mmol scale-up.","evidence_ids":[E[2],T[2]]},
    {"region":"lower-left prose","finding":"Narrative discussion repeats Table 1 results and current-work mechanism summary; matched to existing records instead of duplicated.","evidence_ids":[E[2],T[2]]},
  ],
  3: [
    {"region":"Scheme 5 arrows 3→48, 3→49, 3→50, and 39→51","finding":"Four independent derivatization experiments; 48 is explicitly isolated in one pot, while yield measurement bases for 49-51 are unreported.","evidence_ids":[E[3],T[3]]},
    {"region":"left-column scope discussion","finding":"Narrative summary of products 4-47, the 3 and 39 scale-ups, and entry 14; used as corroborating evidence for existing records.","evidence_ids":[E[3],T[3]]},
    {"region":"right-column asymmetric-version paragraph","finding":"One reported best asymmetric result for product 3 (89% yield, 40% ee); exact catalyst/ligand and measurement methods are referred to unavailable SI and remain unresolved.","evidence_ids":[E[3],T[3]]},
    {"region":"Associated Content and CCDC notice","finding":"SI existence and CCDC 2046574 are contextual notices; SI was not supplied/associated and X-ray characterization is not a separate synthesis.","evidence_ids":[E[3],T[3]]},
  ],
  4: [
    {"region":"Acknowledgments and references 1-8","finding":"No present-paper synthetic experiment, no reaction table, and no characterization peak table; cited literature is bibliography only.","evidence_ids":[E[4],T[4]]},
  ],
  5: [
    {"region":"continuation of references 8-10","finding":"No present-paper synthetic experiment, no reaction scheme, and no characterization peak table; cited allene literature remains bibliography only.","evidence_ids":[E[5],T[5]]},
  ],
}
page_unresolved = {
  1:["Generic Scheme 1 R groups are not concrete identities.","Exact procedures, workups, purifications, and catalyst preparation are in unavailable SI."],
  2:["Dense source-drawn scope structures lack atom-complete verified SMILES.","Scheme 4 loading denominators are not explicitly printed; structured references are marked inferred."],
  3:["Yield measurement bases for 49-51 and the asymmetric 89% result are unreported.","The 40% ee method, major-enantiomer identity, absolute configuration, and exact asymmetric conditions are unavailable.","Transition time for 0 °C to rt in the 39→51 transformation is unreported."],
  4:["Bibliographic citations were not re-extracted as experiments."],
  5:["Bibliographic citations were not re-extracted as experiments."],
}
pages=[]
for p in range(1,6):
    ids=[]
    needles={E[p],T[p]}
    for rx in bundle["reactions"]:
        if needles.intersection(rx["evidence_ids"]):
            ids.append(rx["reaction_id"])
    pages.append({"page":p,"state":"attempted","regions":coverage_regions[p],"record_ids":ids,"unresolved":page_unresolved[p]})
st.write(RUN / "coverage/pages_0001_0005.json", {"document_id":DOC,"semantic_accuracy_verified":False,"pages":pages})

tasks=st.read(RUN / "tasks.json")
for task in tasks:
    if task.get("document_id")==DOC:
        task["state"]="attempted"
        task["reason"]="Page image and extracted text inspected; detailed regions and record mapping are in coverage/pages_0001_0005.json."
st.write(RUN / "tasks.json",tasks)
st.write(RUN / "submitted.json", bundle)
print(json.dumps({"molecules":len(bundle["molecules"]),"reactions":len(bundle["reactions"]),"source_context":len(bundle["source_context"])}, ensure_ascii=False))
