#!/usr/bin/env python3
"""Build the source-specific v3.1 submission for Chem. Commun. 2021, 57, 4783–4786."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_af4daa0432f08df2"
E = {p: [f"{DOC}:page{p}:text", f"{DOC}:page{p}:image"] for p in range(1, 5)}


def f(raw=None, norm=None, ev=(), method="original_extraction", status=None, derivation=None):
    return st.fact(raw, norm, ev, method=method, status=status, derivation=derivation)


def missing(ev=()):
    return st.fact(evidence=ev)


def mol(mid, label, ev, name=None, smiles=None, chirality="unspecified",
        geometry="unspecified", representation=None, assignment=None, mixture=None,
        note=None, axial=None, double_bonds=None):
    m = st.new_molecule(f"{DOC}:{mid}", DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name or label, name or label, ev)
    m["names"]["english"] = f(name or label, name or label, ev)
    m["names"]["abbreviations"] = f(label, [label], ev)
    m["names"]["aliases"] = missing(ev)
    m["cas_numbers"] = missing(ev)
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = f(smiles, smiles, ev, method="structure_derivation",
                         derivation="Candidate graph transcribed from the source drawing/name; parsing is not identity verification.")
        m["generic_template"] = None
        m["stereochemistry"]["isomeric_smiles"] = f(
            None, diag["canonical_isomeric_smiles"], ev, method="structure_derivation",
            derivation="Canonical isomeric SMILES calculated from the unverified candidate graph.")
        m["stereochemistry"]["identity_key"] = f(
            None, diag["stereo_graph_key"], ev, method="structure_derivation",
            derivation="Stereo graph key calculated from the unverified candidate graph.")
        if representation is None:
            representation = "partial"
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(ev)
        m["generic_template"] = None
        m["stereochemistry"]["isomeric_smiles"] = missing(ev)
        m["stereochemistry"]["identity_key"] = missing(ev)
        if representation is None:
            representation = "not_encoded"
    m["stereochemistry"]["atom_centers"] = missing(ev)
    m["stereochemistry"]["double_bonds"] = f(double_bonds, double_bonds, ev) if double_bonds is not None else missing(ev)
    m["stereochemistry"]["axial"] = f(axial, axial, ev) if axial is not None else missing(ev)
    m["stereochemistry"]["mixture"] = f(mixture, mixture, ev) if mixture is not None else missing(ev)
    m["stereochemistry"]["assignment"] = f(assignment, assignment, ev) if assignment is not None else missing(ev)
    m["stereochemistry"]["chirality_status"] = f(chirality, chirality, ev)
    m["stereochemistry"]["geometry_status"] = f(geometry, geometry, ev)
    m["stereochemistry"]["representation_status"] = f(representation, representation, ev)
    m["identity_note"] = f(note, note, ev) if note else missing(ev)
    m["external_properties"] = []
    return m


def template_mol(mid, label, ev, fixed, variables, reason):
    m = mol(mid, label, ev, name=label, chirality="unresolved", geometry="unresolved", representation="not_encoded")
    m["representation_type"] = "generic_template"
    m["smiles"] = missing(ev)
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, ev),
        "variables": [{
            "label": var,
            "attachment_points": missing(ev),
            "definition": f(defn, defn, ev),
            "allowed_values": missing(ev),
            "constraints": missing(ev),
        } for var, defn in variables],
        "instance_mapping": [],
        "reason": f(reason, reason, ev),
    }
    return m


def charge(raw, number, unit, ev, method="original_extraction"):
    return f(raw, {"number": number, "unit": unit}, ev, method=method)


def participant(pid, mid, role, ev, amount=None, equivalents=None, mol_percent=None,
                volume=None, mass=None, ref=None, stage="reaction"):
    p = st.new_participant(pid, f"{DOC}:{mid}", role, ev, stage=stage)
    if amount:
        p["charges"]["amount"] = charge(amount[0], amount[1], "mmol", ev)
    if equivalents:
        p["charges"]["equivalents"] = charge(equivalents[0], equivalents[1], "equiv", ev)
    if mol_percent:
        p["charges"]["mol_percent"] = charge(mol_percent[0], mol_percent[1], "mol%", ev)
    if volume:
        p["charges"]["volume"] = charge(volume[0], volume[1], "mL", ev)
    if mass:
        p["charges"]["mass"] = charge(mass[0], mass[1], mass[2], ev)
    if ref:
        p["equivalence_reference"] = f(ref, f"{DOC}:{ref}", ev)
    return p


def cls(kind, ev):
    if kind == "none":
        return [], [], f("not applicable", "not_applicable", ev), ["No observed transformation; intended bond changes were not assigned."]
    if kind == "unresolved":
        return [], None, f("unresolved", "unresolved", ev), ["Generic family or proposed species lack a safe complete atom mapping."]
    changes = []
    changes.append({"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
                    "before_order": 0, "after_order": 1, "change_type": "formed"})
    changes.append({"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
                    "before_order": 0, "after_order": 1, "change_type": "formed"})
    labels = [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1", "b2"]}]
    if kind == "42":
        changes.append({"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
                        "before_order": 1, "after_order": 2, "change_type": "order_increased"})
        labels.append({"category_number": 13, "code": "DB-C-C", "supporting_change_ids": ["b3"]})
    return ([f(None, c, ev, method="inference", derivation="Class-level net change inferred from the source cycloadduct drawing and reaction description.") for c in changes],
            labels,
            f(None, "partial", ev, method="inference", derivation="Class-level changes only; exact atom identities and full mapping are unresolved."),
            ["Exact atom identities, all bond-order/hydrogen changes, coproducts, and complete atom mapping remain unresolved."])


def reaction(rid, ev, record_type="specific_experiment", categories=None, location=None,
             name=None, transform=None, class_kind="unresolved"):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, ev) if categories is not None else missing(ev)
    r["location"] = f(location, location, ev)
    r["reaction_name"] = f(name, name, ev) if name else missing(ev)
    r["transformation_type"] = f(transform, transform, ev) if transform else missing(ev)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [], "unresolved": []}
    changes, labels, coverage, limits = cls(class_kind, ev)
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = (f(None, labels, ev, method="inference",
        derivation="Labels follow the bond-formation taxonomy from class-level source drawings; no complete atom mapping was asserted.")
        if labels is not None else missing(ev))
    r["classification"]["coverage"] = coverage
    r["classification"]["limitations"] = limits
    r["classification"]["model_votes"] = []
    r["classification"]["mapped_reaction_smiles"] = missing(ev)
    return r


def set_solvent(r, ev, volume=True):
    r["solvent_conditions"]["composition"] = f("THF", [f"{DOC}:thf"], ev)
    if volume:
        r["solvent_conditions"]["total_volume"] = charge("2 mL", 2, "mL", ev)


def add_step(r, ev, temp, hours, participants=None):
    s = st.new_step("reaction", ev)
    s["participant_ids"] = participants or [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    s["operation"] = f("reaction", "reaction", ev)
    s["temperature"] = charge(f"{temp} °C", temp, "°C", ev)
    s["time"] = charge(f"{hours} h", hours, "h", ev)
    r["steps"] = [s]


def outcome(r, oid, kind, targets, raw, norm, ev, basis=None, scope=None, yield_type=None):
    o = st.new_outcome(oid, kind, [f"{DOC}:{x}" for x in targets], ev)
    o["value"] = f(raw, norm, ev)
    if basis:
        o["measurement_basis"] = f(basis, basis, ev)
    if scope:
        o["scope"] = f(scope, scope, ev)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, ev)
    r["outcomes"].append(o)


def common_experiment(r, ev, substrate1="1a", substrate2="2a", product=None, product_kind="32",
                      ligand="DPEPhos", temp=60, catalyst=2.5, ligand_loading=7.5,
                      include_yb=True, include_base=True, degassed=False,
                      amount1=0.4, amount2=0.2, solvent_volume=True, hours=12):
    r["participants"] = [
        participant("substrate_1", substrate1, "reactant", ev, amount=(f"{amount1} mmol", amount1)),
        participant("substrate_2", substrate2, "reactant", ev, amount=(f"{amount2} mmol", amount2)),
        participant("catalyst", "Pd2dba3", "catalyst", ev, mol_percent=(f"{catalyst} mol%", catalyst), ref=substrate2),
        participant("ligand", ligand, "ligand", ev, mol_percent=(f"{ligand_loading} mol%", ligand_loading), ref=substrate2),
    ]
    if include_yb:
        r["participants"].append(participant("lewis_acid", "YbOTf3", "additive", ev, mol_percent=("10 mol%", 10), ref=substrate2))
    if include_base:
        bp = participant("base", "Cs2CO3", "additive", ev, equivalents=("1.0 equiv", 1.0), mol_percent=("100 mol%", 100), ref=substrate2)
        r["participants"].append(bp)
    r["participants"].append(participant("solvent", "thf", "solvent", ev, volume=("2 mL", 2) if solvent_volume else None))
    if product:
        r["participants"].append(participant("product", product, "product", ev))
        r["participants"][-1]["product_role"] = f("desired cycloadduct", "desired", ev)
    set_solvent(r, ev, solvent_volume)
    add_step(r, ev, temp, hours)
    if degassed:
        r["procedure"]["pretreatment"] = f("degassed THF", "degassed THF", ev)
    r["reaction_scale"] = f(f"{amount2} mmol of substrate 2", {"reference_molecule_id": f"{DOC}:{substrate2}", "amount": {"number": amount2, "unit": "mmol"}}, ev)


bundle = json.loads((RUN / "raw_model/standard_init.json").read_text())

# Molecule registry: safe candidates are deliberately sparse; compact products remain source-labelled.
molecules = []
molecules += [
    mol("1a", "1a", E[2], "diethyl 2-vinylidenecyclopropane-1,1-dicarboxylate", "C=C=C1CC1(C(=O)OCC)C(=O)OCC",
        chirality="achiral", geometry="not_applicable", representation="partial",
        axial={"allene": "not stereogenic on candidate local graph because the terminal allene carbon is CH2"},
        note="Candidate terminal-CH2 vinylidenecyclopropane graph; toolkit cumulene support is limited."),
    mol("1b", "1b", E[3], "dimethyl 2-vinylidenecyclopropane-1,1-dicarboxylate", "C=C=C1CC1(C(=O)OC)C(=O)OC",
        chirality="achiral", geometry="not_applicable", representation="partial", axial={"allene": "terminal CH2; no axial stereogenicity on candidate local graph"}),
    mol("1c", "1c", E[3], "diisopropyl 2-vinylidenecyclopropane-1,1-dicarboxylate", "C=C=C1CC1(C(=O)OC(C)C)C(=O)OC(C)C",
        chirality="achiral", geometry="not_applicable", representation="partial", axial={"allene": "terminal CH2; no axial stereogenicity on candidate local graph"}),
    mol("1d", "1d", E[3], "N-propyl-substituted vinylidenecyclopropane-diester 1d", chirality="unspecified", geometry="unspecified",
        axial={"allene": "possible stereogenic allene axis; source gives no assignment"}, note="Compact Scheme 5 drawing not encoded safely."),
    mol("2a", "2a", E[2], "methyl (E)-2-(1-benzyl-2-oxoindolin-3-ylidene)acetate",
        "COC(=O)/C=C1/C(=O)N(Cc2ccccc2)c2ccccc21", chirality="achiral", geometry="specified", representation="encoded",
        double_bonds=[{"descriptor": "E", "basis": "source name and drawing"}], assignment={"alkene": "E"}),
]

sub_names = {
    "2b":"4-bromo methyleneindolinone 2b", "2c":"4-chloro methyleneindolinone 2c",
    "2d":"4-trifluoromethoxy methyleneindolinone 2d", "2e":"4-fluoro methyleneindolinone 2e",
    "2f":"5-bromo methyleneindolinone 2f", "2g":"4,6-dimethyl methyleneindolinone 2g",
    "2h":"N-(4-methylbenzyl) methyleneindolinone 2h", "2i":"N-(3,5-dimethylbenzyl) methyleneindolinone 2i",
    "2j":"N-(2-naphthylmethyl) methyleneindolinone 2j", "2k":"N-(2-thienylmethyl) methyleneindolinone 2k",
    "2l":"N-(methoxycarbonylmethyl) methyleneindolinone 2l", "2m":"N-ethyl methyleneindolinone 2m",
    "2n":"N-isopropyl methyleneindolinone 2n", "2o":"N-(2-methylallyl) methyleneindolinone 2o",
    "2p":"N-cinnamyl methyleneindolinone 2p", "2q":"(E)-1-benzyl-3-(2-oxo-2-phenylethylidene)indolin-2-one 2q",
    "2r":"(E)-1-benzyl-3-(2-oxopropylidene)indolin-2-one 2r",
    "2ph":"N-phenyl methyleneindolinone used for product 4f",
    "2allyl":"N-allyl methyleneindolinone used for product 4j",
}
for sid, name in sub_names.items():
    molecules.append(mol(sid, sid, E[3], name, chirality="achiral", geometry="partially_specified" if sid in ("2q","2r") else "unspecified",
                         representation="not_encoded", double_bonds=[{"descriptor":"E", "basis":"source name"}] if sid in ("2q","2r") else None,
                         note="Compact source structure retained without a forced complete SMILES."))

for prefix, labels in (("3", list("abcdefghijklmnopqrs")), ("4", list("abcdefghijklmn"))):
    for x in labels:
        lab = prefix + x
        if prefix == "3":
            assign = {"relative_stereochemistry": "source-drawn wedge/dash; absolute configuration unassigned"}
            if lab == "3b": assign["structure_assignment"] = "X-ray diffraction; CCDC 2014790 notice, external CIF not supplied"
            molecules.append(mol(lab, lab, E[3] if lab != "3a" else E[2], f"cycloadduct {lab}",
                                 chirality="partially_specified", geometry="unspecified", representation="not_encoded",
                                 assignment=assign, note="Complex spiro product drawing not encoded as a complete stereo SMILES."))
        else:
            assign = {"relative_stereochemistry": "source-drawn major diastereomer; absolute configuration unassigned"}
            if lab == "4b": assign["structure_assignment"] = "X-ray diffraction; CCDC 2024031 notice, external CIF not supplied"
            molecules.append(mol(lab, lab, E[3], f"cycloadduct {lab}", chirality="mixture", geometry="unspecified",
                                 representation="not_encoded", assignment=assign,
                                 mixture={"type":"diastereomeric mixture", "composition":"experiment-specific dr"},
                                 note="Complex spiro product mixture not encoded as a complete stereo SMILES."))

molecules.append(mol("3a_ee31", "3a", E[4], "enantioenriched 3a sample (31% ee)", chirality="mixture", geometry="unspecified",
                     representation="not_encoded", mixture={"type":"enantiomeric mixture", "ee":"31%; major absolute configuration unassigned"},
                     assignment={"absolute_configuration":"unassigned"}))

simple = [
    ("Pd2dba3","Pd2(dba)3"),("YbOTf3","Yb(OTf)3"),("Cs2CO3","Cs2CO3"),("thf","THF"),
    ("XantPhos","XantPhos"),("racBINAP","rac-BINAP"),("XPhos","XPhos"),("DPEPhos","DPEPhos"),
    ("dppe","dppe"),("dppp","dppp"),("dppb","dppb"),("ligand_entry13_conflict","entry-13 ligand: table dppp / prose dppb"),("L1","chiral phosphine ligand L1")]
for mid, name in simple:
    ch = "racemic" if mid == "racBINAP" else "specified" if mid == "L1" else "unspecified"
    molecules.append(mol(mid, name, E[2] if mid != "L1" else E[4], name, chirality=ch, geometry="not_applicable", representation="not_encoded",
                         assignment={"source":"explicit chiral ligand drawing"} if mid == "L1" else None,
                         mixture={"type":"racemic atropisomeric ligand"} if mid == "racBINAP" else None))

molecules += [
    template_mol("prior_propargyl_family", "Scheme 1 propargylic palladium reaction family", E[1],
                 "Propargylic donor A and palladium intermediates B-D with generic R/X/PdL_n/Nu variables",
                 [("R","source-defined substituent"),("X","ester or halide"),("Nu","hard/soft nucleophile")],
                 "Scheme 1 summarizes multiple cited reaction modes rather than one fully specified molecule."),
    template_mol("generic_1", "vinylidenecyclopropane-diester 1", E[2],
                 "Vinylidenecyclopropane bearing gem-diesters and R groups", [("R1","ester alkyl group"),("R2","terminal allene substituent")],
                 "Current-work R-group family in Scheme 2b."),
    template_mol("generic_2", "methyleneindolinone 2", E[2],
                 "Methyleneindolinone bearing R3/R4/R5 variables", [("R3","exocyclic carbonyl substituent"),("R4","benzene substituent"),("R5","nitrogen protecting group")],
                 "Current-work R-group family in Scheme 2b."),
    template_mol("generic_3", "(3+2) cycloadduct family 3", E[2], "Spirooxindole fused to a five-membered ring",
                 [("R1-R5","inherited substrate variables")], "Generic current-work product family."),
    template_mol("generic_4", "(4+2) cycloadduct family 4", E[2], "Spirooxindole fused to a six-membered ring",
                 [("R1-R5","inherited substrate variables")], "Generic current-work product family."),
]
bundle["molecules"] = molecules

reactions = []
prior = reaction("scheme1_prior_modes", E[1], "literature_example", None, "Scheme 1", "different reaction modes of propargyl palladium", "cited propargyl-palladium reaction families", "unresolved")
prior["participants"] = [participant("family", "prior_propargyl_family", "reactant", E[1])]
prior["review"]["unresolved"] = ["Scheme 1 groups cited literature modes without concrete charges or a single experimental arrow."]
reactions.append(prior)

gen = reaction("scheme2_current_generic", E[2], "generic_scheme", None, "Scheme 2b", "palladium-catalyzed divergent cycloadditions", "ligand-switched (3+2)/(4+2) cycloaddition", "unresolved")
gen["participants"] = [participant("vdcp", "generic_1", "reactant", E[2]), participant("indolinone", "generic_2", "reactant", E[2]),
                       participant("p32", "generic_3", "product", E[2]), participant("p42", "generic_4", "product", E[2])]
gen["review"]["unresolved"] = ["Generic R-group family; no single experiment or yield is asserted."]
reactions.append(gen)

table = [
    (1,"XantPhos",80,"3a",11,None,False,True,5,13,False),
    (2,"XantPhos",80,"3a",36,None,True,True,5,13,False),
    (3,"racBINAP",80,None,None,None,True,True,5,13,False),
    (4,"XPhos",80,None,None,None,True,True,5,13,False),
    (5,"DPEPhos",80,"3a",42,None,True,True,5,13,False),
    (6,"dppe",80,"4a",54,"5:1",True,True,5,13,False),
    (7,"dppp",80,"4a",53,"7:1",True,True,5,13,False),
    (8,"DPEPhos",60,"3a",50,None,True,True,5,13,False),
    (9,"DPEPhos",60,"3a",32,None,True,False,5,13,False),
    (10,"DPEPhos",60,"3a",63,None,True,True,2.5,7.5,True),
    (11,"dppe",60,"4a",55,"4:1",True,True,5,13,False),
    (12,"dppp",60,"4a",52,"8:1",True,True,5,13,False),
    (13,"ligand_entry13_conflict",60,"4a",47,"9:1",True,True,5,13,False),
    (14,"dppp",80,"4a",56,"11:1",True,True,2.5,7.5,False),
]
for entry, ligand, temp, product, yld, dr, base, yb, cat, ligload, degas in table:
    ev = E[2] + (E[3] if entry in (10,14) else [])
    r = reaction(f"table1_entry{entry}", ev, categories=["optimization"] + (["substrate_scope"] if entry in (10,14) else []),
                 location=f"Table 1 entry {entry}" + ("; Scheme 3 product 3a" if entry==10 else "; Scheme 4 product 4a" if entry==14 else ""),
                 name="Pd-catalyzed divergent cycloaddition", transform="(3+2) cycloaddition" if product=="3a" else "(4+2) cycloaddition" if product=="4a" else "optimization outcome",
                 class_kind="32" if product=="3a" else "42" if product=="4a" else "none")
    common_experiment(r, ev, product=product, product_kind="32" if product=="3a" else "42", ligand=ligand, temp=temp,
                      catalyst=cat, ligand_loading=ligload, include_yb=yb, include_base=base, degassed=degas)
    r["general_scheme_ids"] = [f"{DOC}:scheme2_current_generic"]
    if yld is not None:
        outcome(r, "yield", "yield", [product], f"{yld}%", {"number":yld,"unit":"%"}, ev, basis="isolated yield", yield_type="isolated")
    if dr:
        outcome(r, "dr", "selectivity", [product], dr, dr, ev, scope={"metric":"dr", "order":"reported major:minor; configurations unassigned"})
    if entry == 1:
        outcome(r, "single_diastereomer", "selectivity", ["3a"], "single diastereomer", "single diastereomer", E[2], scope={"metric":"diastereoselectivity", "configuration":"unassigned"})
    if entry == 3:
        outcome(r, "complex", "other", [], "Complex", "complex product mixture", ev)
    if entry == 4:
        outcome(r, "nr", "no_reaction", [], "NR", "no reaction", ev)
    if entry == 13:
        conflict = f(["table: dppp", "prose: dppb"], None, E[2], status="conflict")
        conflict["alternatives"] = [f("dppp", "dppp", E[2]), f("dppe, dppp and dppb for entries 11–13", "dppb", E[2])]
        r["notes"].append(conflict)
        r["review"]["unresolved"].append("Entry 13 ligand conflicts: table prints dppp; prose identifies entries 11–13 as dppe, dppp and dppb.")
    if entry in (10,14):
        r["notes"].append(f("Same experiment repeated in the corresponding substrate-scope scheme; evidence merged, not double-counted.",
                            "cross-region duplicate merged", ev, method="inference", derivation="Reactants, conditions, product, yield and dr match."))
    reactions.append(r)

scope3 = {
    "3b":("2b",62),"3c":("2c",56),"3d":("2d",57),"3e":("2e",57),"3f":("2f",52),"3g":("2g",58),
    "3h":("2h",57),"3i":("2i",63),"3j":("2j",58),"3k":("2k",56),"3l":("2l",48),"3m":("2m",42),
    "3n":("2n",45),"3o":("2o",56),"3p":("2p",53),"3q":("2q",46),
}
for product,(sub,yld) in scope3.items():
    r = reaction(f"scheme3_{product}", E[3], categories=["substrate_scope"], location=f"Scheme 3 product {product}",
                 name="(3+2) cycloaddition", transform="(3+2) cycloaddition", class_kind="32")
    common_experiment(r,E[3],substrate2=sub,product=product,ligand="DPEPhos",temp=60,catalyst=2.5,ligand_loading=7.5,degassed=True)
    r["general_scheme_ids"]=[f"{DOC}:scheme2_current_generic"]
    outcome(r,"yield","yield",[product],f"{yld}%",{"number":yld,"unit":"%"},E[3],basis="isolated yield",yield_type="isolated")
    reactions.append(r)
r = reaction("scheme3_2r_complex", E[3], categories=["substrate_scope"], location="Scheme 3 substrate 2r",
             name="attempted (3+2) cycloaddition", transform="nonproductive scope attempt", class_kind="none")
common_experiment(r,E[3],substrate2="2r",product=None,ligand="DPEPhos",temp=60,catalyst=2.5,ligand_loading=7.5,degassed=True)
outcome(r,"complex","other",[],"complex","complex product mixture",E[3])
r["review"]["unresolved"]=["Desired product identity, conversion, recovery, byproducts and detection method are unreported."]
reactions.append(r)

scope4 = {
    "4b":("2b",55,"9:1"),"4c":("2c",62,"6:1"),"4d":("2f",57,"9:1"),"4e":("2g",48,"6:1"),
    "4f":("2ph",38,"7:1"),"4g":("2k",53,"9:1"),"4h":("2j",52,"8:1"),"4i":("2p",51,"4:1"),
    "4j":("2allyl",44,"7:1"),"4k":("2q",57,"7:1"),
}
for product,(sub,yld,dr) in scope4.items():
    r = reaction(f"scheme4_{product}", E[3], categories=["substrate_scope"], location=f"Scheme 4 product {product}",
                 name="(4+2) cycloaddition", transform="(4+2) cycloaddition", class_kind="42")
    common_experiment(r,E[3],substrate2=sub,product=product,ligand="dppp",temp=80,catalyst=2.5,ligand_loading=7.5,degassed=True)
    r["general_scheme_ids"]=[f"{DOC}:scheme2_current_generic"]
    outcome(r,"yield","yield",[product],f"{yld}%",{"number":yld,"unit":"%"},E[3],basis="isolated yield",yield_type="isolated")
    outcome(r,"dr","selectivity",[product],dr,dr,E[3],scope={"metric":"dr","order":"reported major:minor; configurations unassigned"})
    if product in ("4f", "4j"):
        r["review"]["unresolved"].append("The compact scheme defines this N-protecting group through the product/source drawing; a standalone substrate label is not printed.")
    reactions.append(r)
r = reaction("scheme4_2r_complex", E[3], categories=["substrate_scope"], location="Scheme 4 substrate 2r",
             name="attempted (4+2) cycloaddition", transform="nonproductive scope attempt", class_kind="none")
common_experiment(r,E[3],substrate2="2r",product=None,ligand="dppp",temp=80,catalyst=2.5,ligand_loading=7.5,degassed=True)
outcome(r,"complex","other",[],"complex","complex product mixture",E[3])
r["review"]["unresolved"]=["Desired product identity, conversion, recovery, byproducts and detection method are unreported."]
reactions.append(r)

scope5 = [
    ("3r","1b","DPEPhos",60,54,None,"32"),("3s","1c","DPEPhos",60,52,None,"32"),
    (None,"1d","DPEPhos",60,None,None,"none"),("4l","1b","dppp",80,56,"8:1","42"),
    ("4m","1c","dppp",80,55,"7:1","42"),("4n","1d","dppp",80,47,"3:1","42"),
]
for idx,(product,sub,lig,temp,yld,dr,ck) in enumerate(scope5,1):
    rid = f"scheme5_{product}" if product else "scheme5_1d_conditionA_complex"
    r = reaction(rid,E[3],categories=["substrate_scope"],location=f"Scheme 5 arrow {idx}",name="divergent cycloaddition",
                 transform="(3+2) cycloaddition" if ck=="32" else "(4+2) cycloaddition" if ck=="42" else "nonproductive scope attempt",class_kind=ck)
    common_experiment(r,E[3],substrate1=sub,substrate2="2a",product=product,ligand=lig,temp=temp,catalyst=2.5,ligand_loading=7.5,degassed=True)
    if yld is not None:
        outcome(r,"yield","yield",[product],f"{yld}%",{"number":yld,"unit":"%"},E[3],basis="isolated yield",yield_type="isolated")
    if dr:
        outcome(r,"dr","selectivity",[product],dr,dr,E[3],scope={"metric":"dr","order":"reported major:minor; configurations unassigned"})
    if product is None:
        outcome(r,"complex","other",[],"complex","complex product mixture",E[3])
        r["review"]["unresolved"]=["Condition A gives a complex mixture; desired 3t, conversion, recovery and byproducts are not established."]
    reactions.append(r)

mech = reaction("scheme6_proposed_mechanism", E[3]+E[4], "proposed_mechanism", None, "Scheme 6",
                "proposed mechanism for divergent cycloadditions", "author-proposed Pd π-propargyl mechanism", "unresolved")
mech["participants"]=[participant("one","1a","reactant",E[3]+E[4]),participant("two","2a","reactant",E[3]+E[4]),
                      participant("three","3a","product",E[3]+E[4]),participant("four","4a","product",E[3]+E[4])]
mech["review"]["unresolved"]=["Intermediates I–VI and paths a/b are proposed, not isolated executed experiments; exact atom mapping is unresolved."]
reactions.append(mech)

asym = reaction("scheme7_asymmetric_3a",E[4],categories=["other"],location="Scheme 7a",name="catalytic asymmetric (3+2) cycloaddition",
                transform="asymmetric (3+2) cycloaddition",class_kind="32")
common_experiment(asym,E[4],product="3a_ee31",ligand="L1",temp=60,catalyst=2.5,ligand_loading=7.5,solvent_volume=False,hours=24)
asym["reaction_scale"] = missing(E[4])
outcome(asym,"yield","yield",["3a_ee31"],"33%",{"number":33,"unit":"%"},E[4])
outcome(asym,"ee","selectivity",["3a_ee31"],"31% ee",{"number":31,"unit":"%"},E[4],scope={"metric":"ee","major absolute configuration":"unassigned"})
asym["review"]["unresolved"]=["Substrate scale, solvent volume, yield method, ee method and absolute configuration are unreported on page 4; Table S1 was not supplied."]
reactions.append(asym)

for product,lig,temp,mass,yld,dr in [("3a","DPEPhos",60,0.65,52,None),("4a","dppp",80,0.58,47,"7:1")]:
    r = reaction(f"scheme7_scaleup_{product}",E[4],categories=["scale_up"],location=f"Scheme 7b scale-up to {product}",
                 name="scale-up divergent cycloaddition",transform="(3+2) cycloaddition" if product=="3a" else "(4+2) cycloaddition",class_kind="32" if product=="3a" else "42")
    common_experiment(r,E[4],product=product,ligand=lig,temp=temp,catalyst=2.5,ligand_loading=7.5,amount1=5.0,amount2=2.5,degassed=True)
    r["participants"][-1]["charges"]["mass"] = charge(f"{mass} g",mass,"g",E[4])
    outcome(r,"product_amount","product_amount",[product],f"{mass} g",{"number":mass,"unit":"g"},E[4])
    outcome(r,"yield","yield",[product],f"{yld}%",{"number":yld,"unit":"%"},E[4])
    if dr:
        outcome(r,"dr","selectivity",[product],dr,dr,E[4],scope={"metric":"dr","order":"reported major:minor; configurations unassigned"})
    r["review"]["unresolved"]=["Yield measurement basis is not reported for Scheme 7b; workup and purification are unavailable without SI."]
    reactions.append(r)

bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id":f"{DOC}:paper_metadata","document_id":DOC,"kind":"metadata","evidence_ids":E[1],
     "content":"Title, authors, DOI 10.1039/d1cc01453f, abstract and publication metadata."},
    {"context_id":f"{DOC}:crystallography","document_id":DOC,"kind":"structure_assignment","evidence_ids":E[1]+E[3],
     "content":"3b and 4b are assigned by X-ray; CCDC 2014790 and 2024031 are external notices only. CIF files were not supplied."},
    {"context_id":f"{DOC}:si_notice","document_id":DOC,"kind":"source_limitation","evidence_ids":E[1]+E[4],
     "content":"The paper advertises ESI and refers to Table S1/Schemes S1-S2, but no SI file was provided or associated."},
    {"context_id":f"{DOC}:characterization_inventory","document_id":DOC,"kind":"negative_inventory","evidence_ids":E[1]+E[2]+E[3]+E[4],
     "content":"No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block occurs on pages 1-4."},
]

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2)+"\n")

stereo = {}
for m in molecules:
    smi = st.value(m["smiles"])
    if smi:
        stereo[m["molecule_id"]] = st.stereo_inspect(smi)
(RUN / "raw_model/stereo_checks.json").write_text(json.dumps(stereo,ensure_ascii=False,indent=2)+"\n")

inventory = {
    "document_id":DOC,
    "pages":[
        {"page":1,"inventory":["title/abstract/metadata","Scheme 1 cited propargyl-palladium reaction modes","ESI/CCDC availability notice"],"reaction_ids":[f"{DOC}:scheme1_prior_modes"]},
        {"page":2,"inventory":["Scheme 2 current generic divergent cycloaddition","Table 1 entries 1-14","scope discussion"],"reaction_ids":[f"{DOC}:scheme2_current_generic"]+[f"{DOC}:table1_entry{i}" for i in range(1,15)]},
        {"page":3,"inventory":["Scheme 3: 3a-3q plus 2r complex","Scheme 4: 4a-4k plus 2r complex","Scheme 5: six arrows","Scheme 6 discussion begins"],
         "reaction_ids":[f"{DOC}:table1_entry10"]+[f"{DOC}:scheme3_3{x}" for x in "bcdefghijklmnopq"]+[f"{DOC}:scheme3_2r_complex",f"{DOC}:table1_entry14"]+[f"{DOC}:scheme4_4{x}" for x in "bcdefghijk"]+[f"{DOC}:scheme4_2r_complex"]+[f"{DOC}:scheme5_{x}" for x in ["3r","3s"]]+[f"{DOC}:scheme5_1d_conditionA_complex"]+[f"{DOC}:scheme5_{x}" for x in ["4l","4m","4n"]]+[f"{DOC}:scheme6_proposed_mechanism"]},
        {"page":4,"inventory":["Scheme 6 proposed mechanism","Scheme 7a asymmetric experiment","Scheme 7b two scale-up arrows","summary/references"],
         "reaction_ids":[f"{DOC}:scheme6_proposed_mechanism",f"{DOC}:scheme7_asymmetric_3a",f"{DOC}:scheme7_scaleup_3a",f"{DOC}:scheme7_scaleup_4a"]},
    ],
    "deduplication":["Table 1 entry 10 = Scheme 3 product 3a","Table 1 entry 14 = Scheme 4 product 4a"],
    "conflicts":["Table 1 entry 13 ligand: table dppp; prose sequence says dppb"],
    "limitations":["No SI supplied","Compact structures intentionally left unencoded where unsafe","Yield/ee methods absent for Scheme 7"],
}
(RUN / "raw_model/page_extraction.json").write_text(json.dumps(inventory,ensure_ascii=False,indent=2)+"\n")

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1,
            "state": "attempted",
            "regions": [
                {"region":"title, abstract, introduction and metadata","finding":"Current paper metadata and overview; no concrete current-paper experiment arrow.","evidence_ids":E[1]},
                {"region":"Scheme 1","finding":"Cited propargyl-palladium reaction modes retained as a literature-example family, not current experiments.","evidence_ids":E[1]},
                {"region":"ESI/CCDC notice","finding":"ESI and CCDC 2014790/2024031 are advertised, but no associated SI/CIF was supplied.","evidence_ids":E[1]},
            ],
            "record_ids": inventory["pages"][0]["reaction_ids"],
            "unresolved": ["No SI or crystallographic file was supplied; the cited literature families lack concrete charges on this page."],
        },
        {
            "page": 2,
            "state": "attempted",
            "regions": [
                {"region":"Scheme 2a-b","finding":"Working hypothesis plus the current generic ligand-switched (3+2)/(4+2) transformation.","evidence_ids":E[2]},
                {"region":"Table 1 rows 1-14 and footnotes a-f","finding":"All fourteen optimization rows recorded separately; rows 3/4 preserve complex/NR, and footnote omissions/loading/degassing are structured.","evidence_ids":E[2]},
                {"region":"optimization and scope prose","finding":"Confirms entry-specific outcomes and exposes the entry-13 table/prose ligand conflict; begins Scheme-3 scope discussion.","evidence_ids":E[2]},
            ],
            "record_ids": inventory["pages"][1]["reaction_ids"],
            "unresolved": ["Entry 13 table prints dppp while prose names dppe/dppp/dppb for entries 11-13; no value was selected."],
        },
        {
            "page": 3,
            "state": "attempted",
            "regions": [
                {"region":"Scheme 3 and caption","finding":"All products 3a-3q and the 2r complex-mixture attempt inventoried; 3a merged with Table 1 entry 10.","evidence_ids":E[3]},
                {"region":"Scheme 4 and caption","finding":"All products 4a-4k and the separate 2r complex-mixture attempt inventoried; yields, dr values and isolated-yield footnote retained; 4a merged with Table 1 entry 14.","evidence_ids":E[3]},
                {"region":"Scheme 5","finding":"Six distinct condition-A/B arrows recorded, including the 1d/condition-A complex mixture and products 3r, 3s, 4l, 4m and 4n.","evidence_ids":E[3]},
                {"region":"mechanistic discussion leading to Scheme 6","finding":"Author-proposed mechanism context; no intermediates were treated as executed experiments.","evidence_ids":E[3]},
            ],
            "record_ids": inventory["pages"][2]["reaction_ids"],
            "unresolved": ["Most compact scope structures lack safe complete stereo SMILES; complex-mixture composition and failed-substrate conversion/recovery are unreported."],
        },
        {
            "page": 4,
            "state": "attempted",
            "regions": [
                {"region":"Scheme 6","finding":"One proposed mechanism record covering paths a/b and intermediates I-VI; not counted as executed experiments.","evidence_ids":E[4]},
                {"region":"Scheme 7a","finding":"One asymmetric experiment to enantioenriched 3a (33% yield, 31% ee); Table S1 is unavailable.","evidence_ids":E[4]},
                {"region":"Scheme 7b","finding":"Two separate 5.0/2.5 mmol scale-up experiments to 3a and 4a with product masses, yields, and 4a dr.","evidence_ids":E[4]},
                {"region":"summary, acknowledgements and references","finding":"No additional current-paper experiment or standalone spectral/mass table.","evidence_ids":E[4]},
            ],
            "record_ids": inventory["pages"][3]["reaction_ids"],
            "unresolved": ["Scheme 7 yield and ee measurement methods, detailed workup/purification, and absolute configuration are unreported on the supplied page."],
        },
    ],
}
(RUN / "coverage/pages_0001_0004.json").write_text(json.dumps(coverage,ensure_ascii=False,indent=2)+"\n")

print(json.dumps({"molecules":len(molecules),"reactions":len(reactions)},ensure_ascii=False))
