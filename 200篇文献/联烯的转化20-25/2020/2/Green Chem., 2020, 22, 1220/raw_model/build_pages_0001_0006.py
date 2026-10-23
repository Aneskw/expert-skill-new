#!/usr/bin/env python3
"""Build evidence-linked v3.1 records for Green Chem. 2020, 22, 1220, PDF pages 1-6.

This is a reproducible source-extraction artifact, not an automatic chemistry
certification.  Every output reaction remains partial/needs_review.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SKILL_SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SKILL_SCRIPTS))
import standard_records as st

DOC = "doc_41522308b99275e3"
TITLE = "Thermal [2 + 2]-cycloaddition between silylalkynes and allenylphenols followed by the nucleophilic addition of water: metal-free and economical synthesis of arylcyclobutenals"


def ev(*pages: int) -> list[str]:
    return [f"{DOC}:page{p}:text" for p in pages]


def F(raw=None, normalized=None, pages=(), method="original_extraction", status=None,
      derivation=None, alternatives=()):
    return st.fact(raw, normalized, ev(*pages), method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def missing():
    return st.fact()


def inferred(normalized, pages, derivation, raw=None):
    return F(raw, normalized, pages, method="inference", derivation=derivation)


def derived(raw, normalized, pages, derivation):
    return F(raw, normalized, pages, method="structure_derivation", derivation=derivation)


def conflict(raw_values, pages, alternatives):
    return F(raw_values, None, pages, status="conflict", alternatives=alternatives)


bundle = json.loads((RUN / "submitted.json").read_text())
dataset = json.loads((RUN / "dataset.json").read_text())
dataset["documents"][0]["title"] = TITLE
bundle["documents"] = copy.deepcopy(dataset["documents"])
st.write(RUN / "dataset.json", dataset)
st.write(RUN / "raw_model/init_pages_0001_0006.json", st.empty(bundle["documents"]))

molecules = {}


def base_molecule(mid, label, name, pages, representation="unresolved"):
    m = st.new_molecule(mid, DOC, ev(*pages))
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, pages)
    m["names"]["original"] = F(name, name, pages)
    m["names"]["english"] = F(name, name, pages)
    m["cas_numbers"] = missing()
    m["representation_type"] = representation
    m["identity_note"] = F(
        "Source label/name retained; CAS not reported.",
        "Source identity retained without external database enrichment.", pages)
    molecules[mid] = m
    return m


def concrete(mid, label, name, smiles, pages, chirality="achiral", geometry="not_applicable",
             representation="not_applicable", assignment=None):
    m = base_molecule(mid, label, name, pages, "smiles")
    m["smiles"] = derived(
        {"source": "structure drawing and/or systematic name", "candidate_smiles": smiles},
        smiles, pages,
        "Candidate connectivity transcribed from the source drawing/name; RDKit parsing is a graph diagnostic, not source verification.")
    diag = st.stereo_inspect(smiles)
    m["stereochemistry"]["chirality_status"] = inferred(
        chirality, pages, "Status assigned from the depicted connectivity and absence/presence of explicitly specified stereochemical marks.")
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, pages, "Status assigned from the depicted connectivity and source stereochemical marks.")
    m["stereochemistry"]["representation_status"] = inferred(
        representation, pages, "Candidate SMILES representation was compared with toolkit potential-stereo diagnostics.")
    if assignment is not None:
        m["stereochemistry"]["assignment"] = F(assignment, assignment, pages)
    if chirality == "unspecified":
        m["stereochemistry"]["atom_centers"] = derived(
            "No wedge/dash or absolute descriptor is printed for the depicted tetrahedral center.",
            [x for x in diag["atom_centers"] if x["CIP"] == "?"], pages,
            "RDKit flags potential tetrahedral centers in the candidate graph; the source leaves their configuration unspecified.")
    else:
        m["stereochemistry"]["atom_centers"] = derived(
            "No stereogenic atom center is specified in the source drawing.", diag["atom_centers"], pages,
            "Toolkit graph diagnostic on the candidate structure; absence of a toolkit flag does not certify all stereochemical possibilities.")
    m["stereochemistry"]["double_bonds"] = derived(
        "No E/Z descriptor is printed.", diag["double_bonds"], pages,
        "Toolkit graph diagnostic; cumulated allene bonds and ring double bonds require separate source interpretation.")
    m["stereochemistry"]["axial"] = derived(
        "Allene axis shown without an axial descriptor." if "=C=" in smiles else "No axial descriptor is printed.",
        [] if chirality == "achiral" else None, pages,
        "The O-allenyl substrates have a terminal CH2 end and are not assigned an allene axial configuration; toolkit axial support is limited.")
    m["stereochemistry"]["mixture"] = missing()
    m["stereochemistry"]["isomeric_smiles"] = missing()
    m["stereochemistry"]["identity_key"] = derived(
        {"candidate_smiles": smiles}, diag["stereo_graph_key"], pages,
        "Generated by standard stereo for graph comparison only; not a sample identity or source verification.")
    return m


def unresolved(mid, label, name, pages, chirality="unresolved", geometry="unresolved",
               representation="unresolved", assignment=None):
    m = base_molecule(mid, label, name, pages)
    m["stereochemistry"]["chirality_status"] = inferred(
        chirality, pages, "The source drawing/name does not support a safely encoded complete stereochemical graph.")
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, pages, "The source drawing/name does not support a safely encoded complete stereochemical graph.")
    m["stereochemistry"]["representation_status"] = inferred(
        representation, pages, "Complete stereochemical representation remains unresolved.")
    if assignment:
        m["stereochemistry"]["assignment"] = F(assignment, assignment, pages)
    return m


def generic(mid, label, name, pages, fixed, variables):
    m = base_molecule(mid, label, name, pages, "generic_template")
    m["generic_template"] = {
        "fixed_structure": F(fixed, fixed, pages),
        "variables": [],
        "instance_mapping": [],
        "reason": F("Generic R-group drawing", "Source uses variable substituents; no unique concrete SMILES is asserted.", pages),
    }
    for v in variables:
        m["generic_template"]["variables"].append({
            "label": v,
            "attachment_points": F(f"Source-labelled attachment point {v}", v, pages),
            "definition": F(f"Variable {v} as printed", None, pages, status="unreadable"),
            "allowed_values": missing(),
            "constraints": missing(),
        })
    for k, val in (("chirality_status", "unresolved"), ("geometry_status", "unresolved"),
                   ("representation_status", "not_encoded")):
        m["stereochemistry"][k] = inferred(val, pages, "Generic variable structure cannot be assigned a single stereochemical state.")
    return m


# Simple agents and solvents.
simple = {
    "water": ("H2O", "water", "O"), "d2o": ("D2O", "deuterium oxide", "[2H]O[2H]"),
    "toluene": ("toluene", "toluene", "Cc1ccccc1"), "benzene": ("benzene", "benzene", "c1ccccc1"),
    "p_xylene": ("p-xylene", "p-xylene", "Cc1ccc(C)cc1"), "methanol": ("MeOH", "methanol", "CO"),
    "ethanol": ("EtOH", "ethanol", "CCO"), "dmf": ("DMF", "N,N-dimethylformamide", "CN(C)C=O"),
    "acetonitrile": ("CH3CN", "acetonitrile", "CC#N"), "dioxane": ("1,4-dioxane", "1,4-dioxane", "O1CCOCC1"),
    "thf": ("THF", "tetrahydrofuran", "C1CCOC1"), "tert_butanol": ("t-BuOH", "tert-butanol", "CC(C)(C)O"),
    "ethyl_acetate": ("AcOEt", "ethyl acetate", "CCOC(C)=O"), "brine": ("brine", "brine", None),
    "sodium_sulfate": ("Na2SO4", "sodium sulfate", "[Na+].[Na+].[O-]S(=O)(=O)[O-]"),
    "ammonium_chloride": ("NH4Cl", "ammonium chloride", "[NH4+].[Cl-]"),
    "diisopropylamine": ("i-Pr2NH", "diisopropylamine", "CC(C)NC(C)C"),
    "cui": ("CuI", "copper(I) iodide", "[Cu]I"),
    "pd_catalyst": ("PdCl2(PPh3)2", "bis(triphenylphosphine)palladium(II) dichloride", None),
    "propargyl_bromide": ("propargyl bromide", "propargyl bromide", "C#CCBr"),
    "potassium_carbonate": ("K2CO3", "potassium carbonate", "[K+].[K+].[O-]C([O-])=O"),
    "potassium_tert_butoxide": ("t-BuOK", "potassium tert-butoxide", "[K+].CC(C)(C)[O-]"),
    "tba_f_trihydrate": ("n-Bu4NF·3H2O", "tetrabutylammonium fluoride trihydrate", None),
    "methyltriphenylphosphonium_bromide": ("Ph3PMeBr", "methyltriphenylphosphonium bromide", None),
    "sodium_chlorite": ("NaClO2", "sodium chlorite", "[Na+].[O-]Cl=O"),
    "sodium_dihydrogen_phosphate": ("NaH2PO4", "sodium dihydrogen phosphate", "[Na+].OP(=O)(O)[O-]"),
    "two_methyl_two_butene": ("2-methyl-2-butene", "2-methyl-2-butene", "CC=C(C)C"),
    "trimethoxybenzene": ("1,3,5-trimethoxybenzene", "1,3,5-trimethoxybenzene", "COc1cc(OC)cc(OC)c1"),
    "n_butyllithium": ("n-BuLi", "n-butyllithium", "[Li]CCCC"),
    "triphenylsilyl_chloride": ("triphenylsilyl chloride", "triphenylsilyl chloride", "Cl[Si](c1ccccc1)(c1ccccc1)c1ccccc1"),
}
for mid, (label, name, smi) in simple.items():
    if smi:
        concrete(mid, label, name, smi, (2, 3, 4, 5, 6), "achiral", "not_applicable", "not_applicable")
    else:
        unresolved(mid, label, name, (3, 4, 5, 6), "unresolved", "unresolved", "unresolved")
unresolved("cu_tc", "CuTc", "copper(I) thiophene-2-carboxylate", (1,), "unresolved", "not_applicable", "unresolved")
unresolved("alumina", "alumina oxide", "alumina", (1,), "not_applicable", "not_applicable", "unresolved")

# Source phenols and alkynes used in general procedure A.
concrete("two_iodophenol", "2-iodophenol", "2-iodophenol", "Oc1ccccc1I", (4, 5))
concrete("two_iodo_four_methylphenol", "2-iodo-4-methylphenol", "2-iodo-4-methylphenol", "Cc1ccc(O)c(I)c1", (5,))
concrete("two_iodo_four_methoxyphenol", "2-iodo-4-methoxyphenol", "2-iodo-4-methoxyphenol", "COc1ccc(O)c(I)c1", (5,))
concrete("four_bromo_two_iodophenol", "4-bromo-2-iodophenol", "4-bromo-2-iodophenol", "Oc1ccc(Br)cc1I", (5,))
concrete("one_iodonaphthalen_two_ol", "1-iodonaphthalen-2-ol", "1-iodonaphthalen-2-ol", "Oc1ccc2ccccc2c1I", (5,))

tips = "[Si](C(C)C)(C(C)C)C(C)C"
alkyne_smiles = {
    "alkyne_tips": f"C#C{tips}", "alkyne_tes": "C#C[Si](CC)(CC)CC",
    "alkyne_tbdms": "C#C[Si](C)(C)C(C)(C)C",
    "alkyne_tbdps": "C#C[Si](c1ccccc1)(c1ccccc1)C(C)(C)C",
    "alkyne_tms": "C#C[Si](C)(C)C", "alkyne_tbutyl": "C#CC(C)(C)C",
    "alkyne_phenyl": "C#Cc1ccccc1", "alkyne_mesityl": "C#Cc1c(C)cc(C)cc1C",
}
for mid, smi in alkyne_smiles.items():
    concrete(mid, mid.replace("alkyne_", "") + " acetylene", mid.replace("alkyne_", "") + " acetylene", smi, (4, 5))

# Concrete O-allenyl-2-ethynyl substrates 1a-o used in Tables 1-3 and Procedures A/B.
core1 = "C=C=COc1ccccc1C#C"
smiles_1 = {
    "1a": core1 + tips,
    "1b": core1 + "[Si](CC)(CC)CC",
    "1c": core1 + "[Si](C)(C)C(C)(C)C",
    "1d": core1 + "[Si](c1ccccc1)(c1ccccc1)C(C)(C)C",
    "1e": core1 + "[Si](c1ccccc1)(c1ccccc1)c1ccccc1",
    "1f": "C=C=COc1ccccc1C#C",
    "1g": "C=C=COc1ccccc1C#CCCCC",
    "1h": "C=C=COc1ccccc1C#CC(C)(C)C",
    "1i": "C=C=COc1ccccc1C#Cc1ccccc1",
    "1j": "C=C=COc1ccccc1C#Cc1c(C)cc(C)cc1C",
    "1k": "C=C=COc1ccccc1C#CI",
    "1l": "C=C=COc1cc(C)ccc1C#C" + tips,
    "1m": "C=C=COc1cc(OC)ccc1C#C" + tips,
    "1n": "C=C=COc1cc(Br)ccc1C#C" + tips,
    "1o": "C=C=COc1ccc2ccccc2c1C#C" + tips,
}
names_1 = {
    "1a": "O-Allenyl-2-[(triisopropylsilyl)ethynyl]phenol",
    "1b": "O-Allenyl-2-[(triethylsilyl)ethynyl]phenol",
    "1c": "O-Allenyl-2-[(tert-butyldimethylsilyl)ethynyl]phenol",
    "1d": "O-Allenyl-2-[(tert-butyldiphenylsilyl)ethynyl]phenol",
    "1e": "O-Allenyl-2-[(triphenylsilyl)ethynyl]phenol",
    "1f": "O-Allenyl-2-ethynylphenol", "1g": "O-Allenyl-2-(n-butylethynyl)phenol",
    "1h": "O-Allenyl-2-(tert-butylethynyl)phenol", "1i": "O-Allenyl-2-(phenylethynyl)phenol",
    "1j": "O-Allenyl-2-(mesitylethynyl)phenol", "1k": "O-Allenyl-2-(iodoethynyl)phenol",
    "1l": "O-Allenyl-4-methyl-2-[(triisopropylsilyl)ethynyl]phenol",
    "1m": "O-Allenyl-4-methoxy-2-[(triisopropylsilyl)ethynyl]phenol",
    "1n": "O-Allenyl-4-bromo-2-[(triisopropylsilyl)ethynyl]phenol",
    "1o": "O-Allenyl-1-[(triisopropylsilyl)ethynyl]naphthalen-2-ol",
}
for label, smi in smiles_1.items():
    pages = (2, 3, 4, 5, 6) if label == "1a" else ((2, 3, 5, 6) if label in {"1f", "1h", "1i", "1j", "1l", "1m", "1n", "1o"} else (2, 4, 5, 6))
    concrete(label, label, names_1[label], smi, pages, "achiral", "not_applicable", "not_applicable")
    molecules[label]["stereochemistry"]["axial"] = derived(
        "No axial descriptor; terminal allenyl CH2 is depicted.", [], pages,
        "The terminal CH2 end has two identical substituents, so this depicted allene is not assigned axial chirality.")

# Cyclobutenal products.  The substituted C4 center is potentially stereogenic but unassigned.
prod_core = "Oc1ccccc1C2=C(C=O)CC2"
smiles_2 = {
    "2a": prod_core + tips, "2b": prod_core + "[Si](CC)(CC)CC",
    "2c": prod_core + "[Si](C)(C)C(C)(C)C",
    "2d": prod_core + "[Si](c1ccccc1)(c1ccccc1)C(C)(C)C",
    "2e": prod_core + "[Si](c1ccccc1)(c1ccccc1)c1ccccc1",
    "2f": "Oc1ccccc1C2=C(C=O)CC2",
    "2h": prod_core + "C(C)(C)C", "2i": prod_core + "c1ccccc1",
    "2l": "Cc1ccc(O)c(C2=C(C=O)CC2" + tips + ")c1",
    "2m": "COc1ccc(O)c(C2=C(C=O)CC2" + tips + ")c1",
    "2n": "Oc1ccc(Br)cc1C2=C(C=O)CC2" + tips,
    "2o": "Oc1ccc2ccccc2c1C3=C(C=O)CC3" + tips,
}
names_2 = {
    "2a": "2-(2-Formyl-4-triisopropylsilylcyclobut-1-enyl)phenol",
    "2b": "2-(2-Formyl-4-triethylsilylcyclobut-1-enyl)phenol",
    "2c": "2-[2-Formyl-4-(tert-butyldimethylsilyl)cyclobut-1-enyl]phenol",
    "2d": "2-[2-Formyl-4-(tert-butyldiphenylsilyl)cyclobut-1-enyl]phenol",
    "2e": "2-(2-Formyl-4-triphenylsilylcyclobut-1-enyl)phenol",
    "2f": "2-(2-Formylcyclobut-1-enyl)phenol", "2h": "2-(2-Formyl-4-tert-butylcyclobut-1-enyl)phenol",
    "2i": "2-(2-Formyl-4-phenylcyclobut-1-enyl)phenol",
    "2l": "2-(2-Formyl-4-triisopropylsilylcyclobut-1-enyl)-4-methylphenol",
    "2m": "2-(2-Formyl-4-triisopropylsilylcyclobut-1-enyl)-4-methoxyphenol",
    "2n": "4-Bromo-2-(2-formyl-4-triisopropylsilylcyclobut-1-enyl)phenol",
    "2o": "1-(2-Formyl-4-triisopropylsilylcyclobut-1-enyl)naphthalene-2-ol",
}
for label, smi in smiles_2.items():
    ch = "achiral" if label == "2f" else "unspecified"
    rep = "not_applicable" if label == "2f" else "partial"
    concrete(label, label, names_2[label], smi, (2, 3, 6), ch, "not_applicable", rep,
             "The source drawings show no wedge/dash or absolute descriptor for the cyclobutene-ring tetrahedral carbon.")

concrete("3a", "3a", "1-Triisopropylsilyl-2H-cyclobuta[c]chromene", "C1Oc2ccccc2C2=C1CC2" + tips, (3, 4), "unspecified", "not_applicable", "partial")
concrete("3j", "3j", "1-Mesityl-2H-cyclobuta[c]chromene", "Cc1cc(C)c(c(C)c1)C1=C2COc3ccccc3C2C1", (2, 3, 6), "unspecified", "not_applicable", "partial")
concrete("4", "4", "2-[4-(Triisopropylsilyl)-2-vinylcyclobut-1-enyl]phenol", "Oc1ccccc1C2=C(C=C)CC2" + tips, (3,), "unspecified", "not_applicable", "partial")
concrete("5", "5", "1-(Triisopropylsilyl)-1,2-dihydro-3H-cyclobuta[c]chromen-3-one", "O=C1Oc2ccccc2C2=C1CC2" + tips, (3,), "unspecified", "not_applicable", "partial")
unresolved("2a_d_labeled", "deuterium-labeled 2a", "deuterium-labeled 2a", (3,),
           "unspecified", "not_applicable", "not_encoded",
           "Scheme 5 depicts D on the substituted cyclobutane carbon and OD at the phenolic oxygen; no stereochemical wedge is shown.")
unresolved("intermediate_II", "II", "zwitterionic/hemiacetal-path intermediate II", (1, 2, 3, 4),
           "unresolved", "unresolved", "not_encoded")

# Generic and prior-work templates.
generic("scheme1_prior_substrate", "Tang substrate", "generic aryl cyclopropane-alkyne substrate", (1,),
        "phenyl-substituted cyclopropane bearing alkynyl-Aryl", ["Aryl"])
generic("scheme1_prior_product", "Tang product", "generic arylcyclobutenal", (1,),
        "phenyl-substituted arylcyclobutenal", ["Aryl"])
generic("scheme2a_substrate", "Scheme 2A substrate", "generic tethered alkyne-allene", (2,),
        "intramolecular alkyne-allene tether through X", ["R1", "R2", "R3", "R4", "X"])
generic("scheme2a_product", "Scheme 2A product", "generic fused cyclobutene", (2,),
        "bi-/tricyclic cyclobutene from intramolecular [2+2]", ["R1", "R2", "R3", "R4", "X"])
generic("scheme2b_substrate", "Jiang substrate", "allenone-ester and aryl alkyne pair", (2,),
        "allenone ester plus aryl alkyne", ["Aryl"])
generic("scheme2b_product", "Jiang product", "polycyclic aromatic product", (2,),
        "water-addition product following acid-promoted [2+2]", ["Aryl"])
generic("generic_1", "1", "generic O-allenyl-2-alkynyl phenol", (2, 4, 6),
        "O-allenyl arene bearing an ortho alkynyl substituent", ["R", "R_prime"])
generic("generic_2", "2", "generic functionalized arylcyclobutenal", (2, 6),
        "arylcyclobutenal with substituents inherited from generic 1", ["R", "R_prime"])
generic("gpA_start", "2-iodophenol derivative", "generic substituted 2-iodophenol", (4,),
        "2-iodophenol framework", ["R"])
generic("gpA_alkyne", "mono-substituted acetylene", "generic mono-substituted acetylene", (4,),
        "terminal acetylene bearing R_prime", ["R_prime"])

reactions = []


def rx(suffix, pages, record_type, categories, name, transformation, location):
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, ev(*pages))
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories, categories, pages) if categories is not None else missing()
    r["location"] = F(location, {"pages": list(pages), "location": location}, pages)
    r["reaction_name"] = F(name, name, pages) if name else missing()
    r["transformation_type"] = F(transformation, transformation, pages) if transformation else missing()
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    reactions.append(r)
    return r


def part(r, pid, mid, role, pages, stage="reaction", mass=None, amount=None, volume=None,
         equiv=None, molpct=None, reference=None, product_role=None, inherited=False):
    p = st.new_participant(pid, mid, role, ev(*pages), stage)
    method = "inherited_general_condition" if inherited else "original_extraction"
    def raw_quantity(value):
        return value.get("raw_text", value) if isinstance(value, dict) else value
    if mass is not None:
        p["charges"]["mass"] = mass if isinstance(mass, dict) and "normalized_value" in mass else F(raw_quantity(mass), mass, pages, method=method)
    if amount is not None:
        p["charges"]["amount"] = F(raw_quantity(amount), amount, pages, method=method)
    if volume is not None:
        p["charges"]["volume"] = F(raw_quantity(volume), volume, pages, method=method)
    if equiv is not None:
        p["charges"]["equivalents"] = F(raw_quantity(equiv), equiv, pages, method=method)
        p["equivalence_reference"] = F(reference, reference, pages, method=method)
    if molpct is not None:
        p["charges"]["mol_percent"] = F(raw_quantity(molpct), molpct, pages, method=method)
        p["equivalence_reference"] = F(reference, reference, pages, method=method)
    if product_role:
        p["product_role"] = F(product_role, product_role, pages)
    r["participants"].append(p)
    return p


def step(r, sid, pages, stage="reaction", operation=None, temperature=None, time=None,
         atmosphere=None, pids=None, inherited=False):
    s = st.new_step(sid, ev(*pages), stage)
    method = "inherited_general_condition" if inherited else "original_extraction"
    s["participant_ids"] = pids or [p["participant_id"] for p in r["participants"] if p["stage"] == stage]
    if operation is not None: s["operation"] = F(operation, operation, pages, method=method)
    if temperature is not None: s["temperature"] = F(temperature.get("raw_text", temperature) if isinstance(temperature, dict) else temperature, temperature, pages, method=method)
    if time is not None: s["time"] = F(time.get("raw_text", time) if isinstance(time, dict) else time, time, pages, method=method)
    if atmosphere is not None: s["atmosphere"] = F(atmosphere, atmosphere, pages, method=method)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, pages, value=None, yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev(*pages))
    if value is not None: o["value"] = F(value.get("raw_text", value) if isinstance(value, dict) else value, value, pages)
    if yield_type is not None: o["yield_type"] = F(yield_type, yield_type, pages)
    if basis is not None: o["measurement_basis"] = F(basis, basis, pages)
    if scope is not None: o["scope"] = F(scope, scope, pages)
    if targets:
        o["product_role"] = inferred("reported target" if kind != "recovery" else "starting material", pages,
                                     "Role follows the source table/product label and participant role.")
    r["outcomes"].append(o)
    return o


def classification(r, pages, changes, limitations):
    if changes == "not_applicable":
        r["classification"]["bond_changes"] = []
        r["classification"]["labels"] = inferred([], pages, "No observed transformation is classified for this no-reaction control.")
        r["classification"]["coverage"] = inferred("not_applicable", pages, "No reaction was reported.")
        r["classification"]["limitations"] = limitations
        return
    if changes is None:
        r["classification"]["bond_changes"] = []
        r["classification"]["labels"] = missing()
        r["classification"]["coverage"] = inferred("unresolved", pages, "Observed decomposition/products are not structurally assigned.")
        r["classification"]["limitations"] = limitations
        return
    facts = []
    labels = {}
    for i, (a, b, before, after, kind) in enumerate(changes, 1):
        ch = {"change_id": f"b{i}", "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
              "before_order": before, "after_order": after, "change_type": kind}
        facts.append(derived("Net bond change inferred from depicted reactant/product structures", ch, pages,
                             "Partial element-level difference; atom correspondence is not claimed complete."))
        if kind in {"formed", "order_increased", "hydrogen_gain"}:
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(f"b{i}")
    r["classification"]["bond_changes"] = facts
    r["classification"]["labels"] = derived(
        "Labels inferred from the depicted net transformation", list(labels.values()), pages,
        "Taxonomy labels include all formed/promoted bonds in the deliberately partial element-level change list.")
    r["classification"]["coverage"] = inferred("partial", pages, "No complete atom mapping is asserted.")
    r["classification"]["limitations"] = limitations


CC = [("C", "C", 0, 1, "formed"), ("C", "C", 0, 1, "formed")]
FULL = CC + [("C", "O", 0, 2, "formed"), ("C", "H", 0, 1, "hydrogen_gain"), ("O", "H", 0, 1, "hydrogen_gain")]


def set_solvent(r, pages, composition, concentration=None, basis=None, inherited=False):
    method = "inherited_general_condition" if inherited else "original_extraction"
    r["solvent_conditions"]["composition"] = F(composition, composition, pages, method=method)
    if concentration is not None: r["solvent_conditions"]["concentration"] = F(concentration, concentration, pages, method=method)
    if basis is not None: r["solvent_conditions"]["concentration_basis"] = F(basis, basis, pages, method=method)


# Page 1 / Scheme 2 literature context and current generic transformation.
r = rx("scheme1_tang_2013", (1,), "literature_example", ["other"],
       "Previously reported synthesis of arylcyclobutenals", "Cu/Ag-catalyzed ring expansion", "Scheme 1")
part(r, "prior_s", "scheme1_prior_substrate", "reactant", (1,)); part(r, "prior_p", "scheme1_prior_product", "product", (1,), product_role="reported prior-work product")
part(r, "cutc", "cu_tc", "catalyst", (1,), molpct={"number": 10, "unit": "mol%"}, reference="scheme1_prior_substrate")
unresolved("ag_otf", "AgOTf", "silver triflate", (1,), "achiral", "not_applicable", "not_applicable")
part(r, "agotf", "ag_otf", "catalyst", (1,), molpct={"number": 10, "unit": "mol%"}, reference="scheme1_prior_substrate")
part(r, "alumina", "alumina", "purification_agent", (1,), stage="purification")
step(r, "metal_catalyzed_rearrangement", (1,), operation="CuTc/AgOTf-catalyzed ring expansion")
step(r, "alumina_treatment", (1,), "purification", operation="alumina oxide treatment", pids=["alumina"])
classification(r, (1,), [("C", "C", 0, 1, "formed")], ["Prior-work generic drawing; substituent identities and complete atom mapping are unresolved."])
r["review"]["unresolved"] = ["Prior-work literature example only; experimental scale, solvent, temperature, time and yield are not reported on this page."]

r = rx("scheme2a_previous_generic", (2,), "literature_example", ["other"], "Intramolecular [2+2] cycloadditions", "thermal or transition-metal-catalyzed alkyne-allene [2+2] cycloaddition", "Scheme 2A")
part(r, "s", "scheme2a_substrate", "reactant", (2,)); part(r, "p", "scheme2a_product", "product", (2,), product_role="prior-work generic product")
step(r, "thermal_or_metal", (2,), operation="heat (Δ) or catalyst M")
classification(r, (2,), CC, ["Generic prior-work scheme with R/X variables; complete mapping unresolved."])
r["review"]["unresolved"] = ["Scheme summarizes cited prior work, not an experiment executed in this paper."]

r = rx("scheme2b_jiang_2018", (2,), "literature_example", ["other"], "Jiang's work (2018)", "acid-promoted [2+2] cycloaddition followed by acid-promoted nucleophilic addition of water", "Scheme 2B")
part(r, "s", "scheme2b_substrate", "reactant", (2,)); part(r, "water", "water", "reagent", (2,)); part(r, "p", "scheme2b_product", "product", (2,), product_role="prior-work generic product")
part(r, "acoh", "water", "additive", (2,))  # placeholder identity corrected below
r["participants"][-1]["molecule_id"] = "acetic_acid"
concrete("acetic_acid", "AcOH", "acetic acid", "CC(=O)O", (2,))
classification(r, (2,), FULL, ["Prior-work generic scheme; element-level classification is partial and not a complete mechanism."])
r["review"]["unresolved"] = ["Scheme summarizes cited prior work; exact scale, temperature, time and yield are absent."]

r = rx("scheme2c_current_generic", (2,), "generic_scheme", ["other"], "Thermal [2+2] cycloaddition followed by nucleophilic addition of water", "one-pot arylcyclobutenal synthesis", "Scheme 2C")
part(r, "s", "generic_1", "reactant", (2,)); part(r, "water", "water", "reagent", (2,)); part(r, "ii", "intermediate_II", "unknown", (2,)); part(r, "p", "generic_2", "product", (2,), product_role="functionalized cyclobutenal")
step(r, "thermal_cycloaddition", (2,), operation="[2+2], Δ"); step(r, "water_addition", (2,), operation="aromatization-promoted nucleophilic addition of H2O")
classification(r, (2,), FULL, ["Generic R/R' structures and elementary atom mapping are unresolved."])
r["review"]["unresolved"] = ["Generic scheme only; use specific experiments for conditions and yields."]


# Table 1 entry 1 is the same experiment as Table 2/Procedure B preparation 2a; create it below.
table1 = [(2, "benzene", 24, 78), (3, "p_xylene", 24, 49), (4, "methanol", 24, 59),
          (5, "ethanol", 24, 63), (6, "dmf", 3, 33), (7, "acetonitrile", 12, 89), (8, "dioxane", 3, 87)]
for entry, solvent_mid, hours, yld in table1:
    r = rx(f"table1_entry{entry}", (2,), "specific_experiment", ["optimization"],
           "Optimization of the reaction conditions", "thermal [2+2]/water-addition cascade", f"Table 1 entry {entry}")
    part(r, "substrate", "1a", "reactant", (2,)); part(r, "solvent", solvent_mid, "solvent", (2,)); part(r, "product", "2a", "product", (2,), product_role="desired cyclobutene")
    set_solvent(r, (2,), simple[solvent_mid][1])
    step(r, "reflux", (2,), operation="reflux", temperature="reflux", time={"number": hours, "unit": "h"})
    outcome(r, "yield_2a", "yield", ["2a"], (2,), {"number": yld, "unit": "%"}, basis=None, scope=f"Table 1 entry {entry}")
    classification(r, (2,), FULL, ["Net element-level changes inferred from the depicted 1a/2a structures; no complete atom mapping."])
    unresolved_items = ["Yield measurement basis is not stated in Table 1.", "Reaction scale, atmosphere and workup are not stated for this optimization row."]
    if entry in (7, 8):
        r["notes"].append(F("Intermediate II was detected by TLC, then underwent water addition after silica-gel purification.", "Intermediate detected by TLC; product formed during/after silica-gel purification.", (2,)))
        unresolved_items.append("The table yield does not separate material formed before versus during silica-gel purification.")
    r["review"]["unresolved"] = unresolved_items


# General procedure B / Table 2 specific experiments.  2a also supplies Table 1 entry 1.
scope = [
    ("1a", "2a", 92, "49 mg", 49, 50, (2, 6), ["optimization", "substrate_scope"]),
    ("1b", "2b", 83, "38 mg", 38, 43, (2, 6), ["substrate_scope"]),
    ("1c", "2c", 73, "37 mg", 37, 43, (2, 6), ["substrate_scope"]),
    ("1d", "2d", 99, "66 mg", 66, 63, (2, 6), ["substrate_scope"]),
    ("1e", "2e", 45, "31 mg", 31, 66, (2, 6), ["substrate_scope"]),
    ("1h", "2h", 62, "23 mg", 23, 34, (2, 3, 6), ["substrate_scope"]),
    ("1i", "2i", 0, None, None, None, (2,), ["substrate_scope"]),
    ("1j", "3j", 80, "35 mg", 35, 44, (2, 3, 6), ["substrate_scope"]),
    ("1l", "2l", 74, "41 mg", 41, 52, (2, 6), ["substrate_scope"]),
    ("1m", "2m", 70, None, None, 55, (2, 3, 6), ["substrate_scope"]),
    ("1n", "2n", 45, None, None, None, (2, 3, 6), ["substrate_scope"]),
    ("1o", "2o", 76, None, None, None, (2, 3, 6), ["substrate_scope"]),
]
for substrate, product, yld, product_mass_raw, product_mass, input_mass, pages, cats in scope:
    suffix = "table1_entry1_table2_2a" if substrate == "1a" else f"table2_{product}"
    r = rx(suffix, pages, "specific_experiment", cats, "General procedure B", "thermal [2+2]/water-addition cascade" if product != "3j" else "thermal [2+2] cycloaddition to stable cyclobutachromene", f"Table 2 {product}; General procedure B where visible")
    gp_b_applies = substrate != "1i"
    if not gp_b_applies:
        mass_fact = None
    elif input_mass is None:
        mass_fact = F("50 mg (general procedure B)", {"number": 50, "unit": "mg"}, (6,), method="inherited_general_condition")
    elif input_mass == 50:
        mass_fact = F("50 mg", {"number": 50, "unit": "mg"}, (6,))
    else:
        mass_fact = conflict(
            ["50 mg (general procedure B)", f"{input_mass} mg (specific preparation)"], (6,),
            [F("50 mg", {"number": 50, "unit": "mg"}, (6,), method="inherited_general_condition"),
             F(f"{input_mass} mg", {"number": input_mass, "unit": "mg"}, (6,))])
    part(r, "substrate", substrate, "reactant", pages, mass=mass_fact,
         amount={"number": 0.16, "unit": "mmol", "raw_text": "0.16 mmol"} if gp_b_applies else None, inherited=gp_b_applies)
    part(r, "solvent", "toluene", "solvent", pages, inherited=gp_b_applies)
    part(r, "product", product, "product", pages, mass={"number": product_mass, "unit": "mg", "raw_text": product_mass_raw} if product_mass is not None else None,
         product_role="stable cyclobutachromene intermediate" if product == "3j" else "reported target product")
    if gp_b_applies:
        set_solvent(r, (6,), "commercially available, not dry toluene", {"number": 0.03, "unit": "M"}, substrate, inherited=True)
        r["reaction_scale"] = F("0.16 mmol", {"reference_molecule_id": substrate, "amount": {"number": 0.16, "unit": "mmol"}}, (6,), method="inherited_general_condition")
        step(r, "purge_and_seal", (6,), operation="place in test tube, purge with nitrogen and seal", atmosphere="nitrogen", inherited=True)
        step(r, "reflux", (2, 6), operation="reflux", temperature="reflux", time={"number": 2, "unit": "h"}, atmosphere="nitrogen", inherited=True)
        step(r, "concentrate", (6,), "workup", operation="cool to room temperature and concentrate in vacuo", inherited=True)
        step(r, "column", (6,), "purification", operation="neutral silica gel column; hexane/AcOEt 10:1 to 5:1", inherited=True)
    else:
        set_solvent(r, (2,), "toluene")
        step(r, "reflux", (2,), operation="reflux", temperature="reflux", time={"number": 2, "unit": "h"})
    basis = "isolated after column chromatography; product mass reported" if product_mass is not None else None
    outcome(r, f"yield_{product}", "yield", [product], pages, {"number": yld, "unit": "%"}, "reported yield", basis, "Table 2/General procedure B")
    classification(r, pages, CC if product == "3j" else FULL, ["Partial net classification; complete atom mapping and proton/water atom origins are unresolved."])
    if gp_b_applies:
        r["notes"].append(F("General procedure B prints 50 mg for compounds 1a-o; specific preparations print substrate-dependent masses.", "General/specific mass discrepancy retained in the substrate mass fact where applicable.", (6,)))
    r["review"]["unresolved"] = (["General-procedure 50 mg versus specific input mass is retained as a conflict where values differ."] if gp_b_applies else ["Only the Table 2 toluene/reflux/2 h conditions are applied; Procedure B is not inherited for the excluded 2i outcome."]) + ["Complete stereochemical assignment of the newly formed tetrahedral ring carbon is not reported."]
    if basis is None: r["review"]["unresolved"].append("Yield measurement basis is not explicitly stated within the viewed pages for this entry.")

# Table 2 decomposition/no-reaction cases (no invented product participants).
for substrate, finding in (("1f", "decomposed"), ("1g", "decomposed")):
    r = rx(f"table2_{substrate}_decomposition", (2, 3), "specific_experiment", ["substrate_scope"], "General procedure B", "failed thermal cascade", f"Table 2 {substrate}; page 3 discussion")
    part(r, "substrate", substrate, "reactant", (2, 3)); part(r, "solvent", "toluene", "solvent", (2, 3))
    set_solvent(r, (2,), "toluene")
    step(r, "reflux", (2,), operation="reflux", temperature="reflux", time={"number": 2, "unit": "h"})
    outcome(r, "decomposition", "decomposition", [substrate], (3,), finding, scope="Table 2 reports 0% target and prose states substrate decomposed")
    classification(r, (2, 3), None, ["Decomposition products are not assigned; no bond labels are invented."])
    r["review"]["unresolved"] = ["Decomposition products, conversion, recovery and detection method are not reported."]

r = rx("table2_1k_no_reaction", (2, 3), "specific_experiment", ["substrate_scope", "control"], "General procedure B", "no reaction", "Table 2 2k / page 3 discussion")
part(r, "substrate", "1k", "reactant", (2, 3)); part(r, "solvent", "toluene", "solvent", (2, 3))
set_solvent(r, (2,), "toluene")
step(r, "reflux", (2,), operation="reflux", temperature="reflux", time={"number": 2, "unit": "h"})
outcome(r, "no_reaction", "no_reaction", ["1k"], (3,), True, scope="reaction did not proceed; Table 2 shows 0%")
classification(r, (2, 3), "not_applicable", ["No product is invented for the explicitly reported no-reaction case."])
r["review"]["unresolved"] = ["Starting-material recovery, observation method and detection limit are not reported."]

# Generic procedure B scheme.
r = rx("general_procedure_B", (6,), "generic_scheme", ["other"], "General procedure B", "generic thermal cascade", "Page 6 general procedure B")
part(r, "s", "generic_1", "reactant", (6,), mass={"number": 50, "unit": "mg"}, amount={"number": 0.16, "unit": "mmol"})
part(r, "solvent", "toluene", "solvent", (6,)); part(r, "p", "generic_2", "product", (6,), product_role="generic product")
set_solvent(r, (6,), "commercially available, not dry toluene", {"number": 0.03, "unit": "M"}, "generic_1")
step(r, "purge", (6,), operation="purge with nitrogen and seal", atmosphere="nitrogen"); step(r, "reflux", (6,), operation="reflux", temperature="reflux", time={"number": 2, "unit": "h"}, atmosphere="nitrogen")
step(r, "concentrate", (6,), "workup", operation="cool to room temperature and concentrate in vacuo"); step(r, "column", (6,), "purification", operation="neutral silica gel; hexane/AcOEt 10:1 to 5:1")
classification(r, (6,), FULL, ["Generic R/R' scheme; exact atom mapping and branch to 3j are unresolved."])
r["review"]["unresolved"] = ["Generic mass 50 mg conflicts with several printed specific substrate masses; specific records retain both values."]

# Scheme 3 derivatizations.
deriv = [
    ("scheme3_desilylation_2f", "2f", 94, [("tba", "tba_f_trihydrate", "reagent", {"number": 3, "unit": "equiv", "raw_text": "3.0 eq."})], "dmf", {"number": 80, "unit": "°C", "raw_text": "80 °C"}, {"number": 4, "unit": "h", "raw_text": "4 h"}, [("C", "H", 0, 1, "hydrogen_gain")], "desilylation"),
    ("scheme3_wittig_4", "4", 76, [("ylide", "methyltriphenylphosphonium_bromide", "reagent", {"number": 2.2, "unit": "equiv", "raw_text": "2.2 eq."}), ("base", "potassium_tert_butoxide", "reagent", {"number": 2.2, "unit": "equiv", "raw_text": "2.2 eq."})], "thf", "r.t.", {"number": 12, "unit": "h", "raw_text": "12 h"}, [("C", "C", 0, 2, "formed")], "Wittig reaction"),
    ("scheme3_pinnick_5", "5", 74, [("chlorite", "sodium_chlorite", "reagent", {"number": 5, "unit": "equiv", "raw_text": "5.0 eq."}), ("phosphate", "sodium_dihydrogen_phosphate", "reagent", {"number": 5, "unit": "equiv", "raw_text": "5.0 eq."}), ("scavenger", "two_methyl_two_butene", "reagent", {"number": 30, "unit": "equiv", "raw_text": "30 eq."})], "tert_butanol", "r.t.", {"number": 1, "unit": "h", "raw_text": "1 h"}, [("C", "O", 0, 1, "formed")], "Pinnick oxidation followed by intramolecular cyclization"),
]
for suffix, product, yld, reagents, solvent_mid, temp, tm, changes, name in deriv:
    r = rx(suffix, (3,), "specific_experiment", ["derivatization"], name, name, "Scheme 3")
    part(r, "substrate", "2a", "reactant", (3,)); part(r, "solvent", solvent_mid, "solvent", (3,)); part(r, "product", product, "product", (3,), product_role="derivatization product")
    for pid, mid, role, eqv in reagents: part(r, pid, mid, role, (3,), equiv=eqv, reference="2a")
    if suffix == "scheme3_pinnick_5": part(r, "water", "water", "solvent", (3,))
    set_solvent(r, (3,), "H2O/t-BuOH" if suffix == "scheme3_pinnick_5" else simple[solvent_mid][1])
    step(r, "reaction", (3,), operation=name, temperature=temp, time=tm)
    outcome(r, f"yield_{product}", "yield", [product], (3,), {"number": yld, "unit": "%"}, basis=None, scope="Scheme 3")
    classification(r, (3,), changes, ["Partial net functional-group classification from Scheme 3; complete reagent/product atom mapping unresolved."])
    r["review"]["unresolved"] = ["Scheme 3 does not state reaction scale or yield measurement basis on the viewed page."]

# Table 3 dry-solvent mechanistic experiments.
table3 = [
    (1, "toluene", 2, None, 15, 5), (2, "p_xylene", 24, None, 9, 3),
    (3, "benzene", 24, None, 69, 3), (4, "methanol", 24, None, 65, 0),
    (5, "ethanol", 24, None, 73, 7), (6, "acetonitrile", 12, None, 68, 9),
    (7, "dioxane", 12, None, 62, 4), (8, "toluene", 2, 1, 0, 92),
    (9, "toluene", 2, 2, 0, 92),
]
for entry, solvent_mid, hours, water_eq, y3, y2 in table3:
    r = rx(f"table3_entry{entry}", (3,), "specific_experiment", ["optimization", "mechanistic_study"], "Reactions in dry solvent", "thermal [2+2] with/without subsequent water addition", f"Table 3 entry {entry}")
    part(r, "substrate", "1a", "reactant", (3,)); part(r, "solvent", solvent_mid, "solvent", (3,)); part(r, "intermediate", "3a", "product", (3,), product_role="cyclobutachromene intermediate"); part(r, "product", "2a", "product", (3,), product_role="water-addition product")
    if water_eq is not None: part(r, "water", "water", "additive", (3,), equiv={"number": water_eq, "unit": "equiv"}, reference="1a")
    part(r, "internal_standard", "trimethoxybenzene", "reagent", (3,), stage="characterization")
    set_solvent(r, (3,), "dry " + simple[solvent_mid][1])
    step(r, "reflux", (3,), operation="reflux in dry solvent", temperature="reflux", time={"number": hours, "unit": "h"}, atmosphere="nitrogen")
    outcome(r, "yield_3a", "yield", ["3a"], (3,), {"number": y3, "unit": "%"}, "NMR yield", "NMR using 1,3,5-trimethoxybenzene as internal standard", f"Table 3 entry {entry}")
    outcome(r, "yield_2a", "yield", ["2a"], (3,), {"number": y2, "unit": "%"}, "NMR yield", "NMR using 1,3,5-trimethoxybenzene as internal standard", f"Table 3 entry {entry}")
    classification(r, (3,), CC if y2 == 0 else FULL, ["Products are both reported outcomes; classification is the union of observed branches and remains partial."])
    r["review"]["unresolved"] = ["Reaction scale and exact dry-solvent preparation are not reported in Table 3.", "NMR observation time and detection limit are not reported."]

# Scale-up Scheme 4, two independent executions.
for cond, solvent_mid, water_eq, yld, atmosphere in (("A", "toluene", 2, 83, "nitrogen"), ("B", "water", None, 68, "air")):
    r = rx(f"scheme4_scaleup_{cond.lower()}", (3,), "specific_experiment", ["scale_up"], "Scale-up experiment", "thermal [2+2]/water-addition cascade", f"Scheme 4 conditions {cond}")
    part(r, "substrate", "1a", "reactant", (3,), amount={"number": 2, "unit": "mmol"}); part(r, "solvent", solvent_mid, "solvent", (3,)); part(r, "product", "2a", "product", (3,), product_role="scale-up product")
    if water_eq is not None: part(r, "water", "water", "additive", (3,), equiv={"number": water_eq, "unit": "equiv"}, reference="1a")
    r["reaction_scale"] = F("1a (2 mmol)", {"reference_molecule_id": "1a", "amount": {"number": 2, "unit": "mmol"}}, (3,))
    set_solvent(r, (3,), "toluene with H2O (2 eq.)" if cond == "A" else "H2O (solvent)")
    step(r, "reflux", (3,), operation="reflux", temperature="reflux", time={"number": 2, "unit": "h"}, atmosphere=atmosphere)
    outcome(r, "yield_2a", "yield", ["2a"], (3,), {"number": yld, "unit": "%"}, basis=None, scope=f"Scheme 4 conditions {cond}")
    classification(r, (3,), FULL, ["Partial net classification; complete atom mapping unresolved."])
    r["review"]["unresolved"] = ["Scheme 4 does not state yield measurement basis or workup details on the viewed page."]

# D2O labeling experiment.
r = rx("scheme5_d2o_labeling", (3,), "specific_experiment", ["mechanistic_study"], "D2O labeling experiment", "thermal [2+2]/D2O-addition labeling", "Scheme 5")
part(r, "substrate", "1a", "reactant", (3,)); part(r, "d2o", "d2o", "reagent", (3,), equiv={"number": 2, "unit": "equiv"}, reference="1a"); part(r, "solvent", "dioxane", "solvent", (3,)); part(r, "labeled", "2a_d_labeled", "product", (3,), product_role="deuterium-labeled product"); part(r, "unlabeled", "2a", "product", (3,), product_role="unlabeled product")
set_solvent(r, (3,), "dioxane")
step(r, "reflux", (3,), operation="reflux with D2O", temperature="reflux", time={"number": 2, "unit": "h"})
outcome(r, "labeled_unlabeled_ratio", "selectivity", ["2a_d_labeled", "2a"], (3,), "1:1", basis=None, scope="deuterium-labeled 2a : unlabeled 2a")
classification(r, (3,), FULL, ["C-D/O-D are classified under C-H/O-H isotope-parent categories; isotope-complete atom mapping is not asserted."])
r["review"]["unresolved"] = ["No yield, measurement method or isotopologue-resolved analytical basis is stated in Scheme 5.", "The drawn D/OD positions are retained textually because a verified isotope/stereo SMILES was not established."]

# Proposed mechanism Scheme 6.
r = rx("scheme6_proposed_mechanism", (3, 4), "proposed_mechanism", ["mechanistic_study"], "Proposed mechanism", "[2+2] cycloaddition, tautomerization, water attack and hemiacetal cleavage", "Scheme 6")
part(r, "substrate", "1a", "reactant", (4,)); part(r, "intermediate3", "3a", "unknown", (4,)); part(r, "ii", "intermediate_II", "unknown", (4,)); part(r, "water", "water", "reagent", (4,)); part(r, "product", "2a", "product", (4,), product_role="mechanistic product")
for sid, op in (("cycloaddition", "thermal [2+2] cycloaddition"), ("tautomerization", "tautomerization"), ("water_attack", "nucleophilic attack of H2O"), ("cleavage", "cleavage of hemiacetal")):
    step(r, sid, (4,), operation=op)
classification(r, (4,), FULL, ["Author-proposed multi-arrow mechanism; not an executed experiment and not a complete elementary atom map."])
r["review"]["unresolved"] = ["Scheme 6 is proposed, not experimentally proven; elementary proton transfers and atom origins remain unresolved."]


# General procedure A and its specific three-step preparations.
def add_gpA_participants(r, start_mid, alkyne_mid, pages, start_mass, start_amount):
    part(r, "phenol", start_mid, "reactant", pages,
         mass={"number": start_mass, "unit": "mg", "raw_text": f"{start_mass} mg"},
         amount={"number": start_amount, "unit": "mmol", "raw_text": f"{start_amount:.2f} mmol"},
         equiv={"number": 1, "unit": "equiv", "raw_text": "1.0 equiv."}, reference=start_mid)
    inherited_pages = (4,)
    part(r, "alkyne", alkyne_mid, "reactant", inherited_pages, equiv={"number": 1.15, "unit": "equiv", "raw_text": "1.15 equiv."}, reference=start_mid, inherited=True)
    part(r, "base_a1", "diisopropylamine", "reagent", inherited_pages, equiv={"number": 3, "unit": "equiv", "raw_text": "3.0 equiv."}, reference=start_mid, inherited=True)
    part(r, "cu", "cui", "catalyst", inherited_pages, molpct={"number": 2, "unit": "mol%", "raw_text": "2 mol%"}, reference=start_mid, inherited=True)
    part(r, "pd", "pd_catalyst", "catalyst", inherited_pages, molpct={"number": 1, "unit": "mol%", "raw_text": "1 mol%"}, reference=start_mid, inherited=True)
    part(r, "thf_a1", "thf", "solvent", inherited_pages, inherited=True)
    part(r, "propargyl", "propargyl_bromide", "reagent", inherited_pages, equiv={"number": 1.5, "unit": "equiv", "raw_text": "1.5 equiv."}, reference=start_mid, inherited=True)
    part(r, "k2co3", "potassium_carbonate", "reagent", inherited_pages, equiv={"number": 2, "unit": "equiv", "raw_text": "2.0 equiv."}, reference=start_mid, inherited=True)
    part(r, "dmf_a2", "dmf", "solvent", inherited_pages, inherited=True)
    part(r, "tbuok", "potassium_tert_butoxide", "reagent", inherited_pages, equiv={"number": 0.5, "unit": "equiv", "raw_text": "0.5 equiv."}, reference=start_mid, inherited=True)
    part(r, "thf_a3", "thf", "solvent", inherited_pages, inherited=True)
    for pid, mid, role in (("nh4cl", "ammonium_chloride", "workup_agent"), ("water", "water", "workup_agent"), ("acoet", "ethyl_acetate", "workup_agent"), ("brine", "brine", "workup_agent"), ("drying", "sodium_sulfate", "workup_agent")):
        part(r, pid, mid, role, inherited_pages, stage="workup", inherited=True)


def add_gpA_steps(r, pages, a3_temp="ambient temperature"):
    reaction_pids = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    step(r, "A1_coupling", (4,), operation="stir Sonogashira mixture in THF (0.2 M)", temperature="ambient temperature", time={"number": 2, "unit": "h"}, atmosphere="nitrogen", pids=reaction_pids, inherited=True)
    step(r, "A1_workup", (4,), "workup", operation="quench with saturated NH4Cl; extract with AcOEt; wash with brine; dry over Na2SO4; concentrate", pids=[p["participant_id"] for p in r["participants"] if p["stage"] == "workup"], inherited=True)
    step(r, "A1_short_column", (4,), "purification", operation="neutral silica gel, hexane/AcOEt 20:1", inherited=True)
    step(r, "A2_propargylation", (4,), operation="stir obtained 2-substituted-ethynylphenol with propargyl bromide/K2CO3 in DMF (0.5 M)", temperature="ambient temperature", time={"number": 6, "unit": "h"}, atmosphere="nitrogen", pids=reaction_pids, inherited=True)
    step(r, "A2_workup", (4,), "workup", operation="quench with water; extract with AcOEt; wash with brine; dry over Na2SO4; concentrate", inherited=True)
    step(r, "A2_short_column", (4,), "purification", operation="neutral silica gel, hexane/AcOEt 50:1", inherited=True)
    modified = isinstance(a3_temp, dict)
    a3_pages = pages if modified else (4,)
    step(r, "A3_isomerization", a3_pages, operation="stir O-propargyl intermediate with t-BuOK in THF (0.2 M)", temperature=a3_temp, time={"number": 1, "unit": "h"}, atmosphere="nitrogen", pids=reaction_pids, inherited=not modified)
    step(r, "A3_workup", (4,), "workup", operation="quench with water; extract with AcOEt; wash with brine; dry over Na2SO4; concentrate", inherited=True)
    step(r, "A3_column", (4,), "purification", operation="neutral flash silica gel, hexane/AcOEt 150:1", inherited=True)


r = rx("general_procedure_A", (4,), "generic_scheme", ["other"], "General procedure A", "Sonogashira coupling, O-propargylation and propargyl-to-allene isomerization", "Page 4 general procedure A")
part(r, "phenol", "gpA_start", "reactant", (4,), equiv={"number": 1, "unit": "equiv"}, reference="gpA_start")
part(r, "alkyne", "gpA_alkyne", "reactant", (4,), equiv={"number": 1.15, "unit": "equiv"}, reference="gpA_start")
part(r, "product", "generic_1", "product", (4,), product_role="generic O-allenyl-2-alkynyl phenol")
add_gpA_participants(r, "gpA_start", "gpA_alkyne", (4,), 0, 0)
# Remove duplicate generic starting participants created by the helper and keep explicit quantities from the scheme.
seen = set(); r["participants"] = [p for p in r["participants"] if not ((p["molecule_id"], p["role"]) in seen or seen.add((p["molecule_id"], p["role"])))]
add_gpA_steps(r, (4,))
classification(r, (4,), [("C", "C", 0, 1, "formed"), ("C", "O", 0, 1, "formed"), ("C", "C", 1, 2, "order_increased"), ("C", "H", 0, 1, "hydrogen_gain")], ["Union of partial stage-level net changes; intermediates and atom mapping are unresolved."])
r["review"]["unresolved"] = ["Generic R/R' identities and intermediate quantities are unspecified.", "Scheme and prose give procedures but no generic isolated yield."]

prepA = [
    ("1a", "two_iodophenol", "alkyne_tips", 550, 2.50, 586, 75, (4,), "ambient temperature"),
    ("1b", "two_iodophenol", "alkyne_tes", 550, 2.50, 122, 18, (4,), "ambient temperature"),
    ("1c", "two_iodophenol", "alkyne_tbdms", 550, 2.50, 360, 53, (4,), "ambient temperature"),
    ("1d", "two_iodophenol", "alkyne_tbdps", 550, 2.50, 710, 72, (4,), "ambient temperature"),
    ("1f", "two_iodophenol", "alkyne_tms", 550, 2.50, 160, 41, (4, 5), "ambient temperature"),
    ("1h", "two_iodophenol", "alkyne_tbutyl", 550, 2.50, 277, 53, (5,), "ambient temperature"),
    ("1i", "two_iodophenol", "alkyne_phenyl", 550, 2.50, 294, 51, (5,), {"number": 0, "unit": "°C"}),
    ("1j", "two_iodophenol", "alkyne_mesityl", 550, 2.50, 466, 68, (5,), {"number": 0, "unit": "°C"}),
    ("1l", "two_iodo_four_methylphenol", "alkyne_tips", 585, 2.50, 530, 65, (5,), "ambient temperature"),
    ("1m", "two_iodo_four_methoxyphenol", "alkyne_tips", 625, 2.50, 608, 71, (5,), "ambient temperature"),
    ("1n", "four_bromo_two_iodophenol", "alkyne_tips", 747, 2.50, 460, 47, (5,), "ambient temperature"),
    ("1o", "one_iodonaphthalen_two_ol", "alkyne_tips", 675, 2.50, 544, 60, (5,), "ambient temperature"),
]
for product, start_mid, alkyne_mid, start_mass, start_amount, product_mass, yld, pages, a3_temp in prepA:
    r = rx(f"preparation_{product}", pages, "specific_experiment", ["other"], "General procedure A", "three-step synthesis of O-allenyl-2-alkynyl phenol", f"Specific preparation of {product}")
    add_gpA_participants(r, start_mid, alkyne_mid, pages, start_mass, start_amount)
    part(r, "product", product, "product", pages, mass={"number": product_mass, "unit": "mg", "raw_text": f"{product_mass} mg"}, product_role="three-step product")
    r["reaction_scale"] = F(f"{start_mid}: {start_amount} mmol", {"reference_molecule_id": start_mid, "amount": {"number": start_amount, "unit": "mmol"}}, pages)
    r["solvent_conditions"]["composition"] = F("THF (A.1, A.3); DMF (A.2)", ["THF", "DMF"], (4,), method="inherited_general_condition")
    add_gpA_steps(r, pages, a3_temp)
    outcome(r, f"overall_yield_{product}", "yield", [product], pages, {"number": yld, "unit": "%", "raw_text": f"{yld}% (3 steps)"}, "three-step isolated yield", "isolated after column chromatography; product mass reported", f"overall from starting iodophenol to {product}")
    outcome(r, f"amount_{product}", "product_amount", [product], pages, {"number": product_mass, "unit": "mg", "raw_text": f"{product_mass} mg"}, basis="isolated mass", scope=f"three-step product {product}")
    classification(r, pages, [("C", "C", 0, 1, "formed"), ("C", "O", 0, 1, "formed"), ("C", "C", 1, 2, "order_increased"), ("C", "H", 0, 1, "hydrogen_gain")], ["Overall three-stage union; no intermediate atom mapping or quantitative intermediate recovery is asserted."])
    r["general_scheme_ids"] = [f"{DOC}:general_procedure_A"]
    r["review"]["unresolved"] = ["Intermediate masses/yields are unreported; the printed yield is retained as a three-step overall yield.", "Cooling/warming transition times and workup volumes are unreported."]
    if product == "1f": r["review"]["unresolved"].append("The source states trimethylsilylacetylene was used in A.1-A.2 and terminal alkyne 1f resulted after A.3; the desilylation detail is not separately quantified.")

# Standalone preparation of 1e from 1f.
r = rx("preparation_1e_from_1f", (5,), "specific_experiment", ["other"], "Preparation of 1e from 1f", "alkynyl C-silylation", "Page 5 preparation of 1e")
part(r, "substrate", "1f", "reactant", (5,), mass={"number": 200, "unit": "mg"}, amount={"number": 1.28, "unit": "mmol"}); part(r, "solvent", "thf", "solvent", (5,)); part(r, "n_buli", "n_butyllithium", "reagent", (5,), equiv={"number": 1, "unit": "equiv"}, reference="1f"); part(r, "silyl_chloride", "triphenylsilyl_chloride", "reagent", (5,), equiv={"number": 2, "unit": "equiv"}, reference="1f"); part(r, "product", "1e", "product", (5,), mass={"number": 397, "unit": "mg"}, product_role="silylated alkyne")
for pid, mid, role in (("nh4cl", "ammonium_chloride", "workup_agent"), ("acoet", "ethyl_acetate", "workup_agent"), ("brine", "brine", "workup_agent"), ("drying", "sodium_sulfate", "workup_agent")):
    part(r, pid, mid, role, (5,), stage="workup")
r["reaction_scale"] = F("1f 1.28 mmol", {"reference_molecule_id": "1f", "amount": {"number": 1.28, "unit": "mmol"}}, (5,))
set_solvent(r, (5,), "THF", {"number": 0.2, "unit": "M"}, "1f")
step(r, "n_buli_addition", (5,), operation="slowly add n-BuLi to 1f in a flame-dried flask", temperature={"number": -78, "unit": "°C"}, atmosphere="nitrogen")
step(r, "lithiation_hold", (5,), operation="stir", temperature={"number": -78, "unit": "°C"}, time={"number": 15, "unit": "min"}, atmosphere="nitrogen")
step(r, "silyl_chloride_addition", (5,), operation="add triphenylsilyl chloride", temperature={"number": -78, "unit": "°C"}, atmosphere="nitrogen")
step(r, "warm_to_room_temperature", (5,), operation="allow the mixture to warm to room temperature", temperature="room temperature", atmosphere="nitrogen")
step(r, "post_addition_hold", (5,), operation="continue reaction", temperature="room temperature", time={"number": 3, "unit": "h"}, atmosphere="nitrogen")
step(r, "workup", (5,), "workup", operation="quench with saturated NH4Cl; extract with AcOEt; wash with brine; dry over Na2SO4; concentrate")
step(r, "column", (5,), "purification", operation="neutral silica gel, hexane/AcOEt 100:1")
outcome(r, "yield_1e", "yield", ["1e"], (5,), {"number": 75, "unit": "%"}, "isolated yield", "isolated after column chromatography; product mass reported", "1f to 1e")
outcome(r, "amount_1e", "product_amount", ["1e"], (5,), {"number": 397, "unit": "mg"}, basis="isolated mass", scope="1e")
classification(r, (5,), [("C", "Si", 0, 1, "formed")], ["Partial net C-Si formation; lithium/chloride coproduct accounting unresolved."])
r["review"]["unresolved"] = ["Warming-ramp duration is not reported; it is not merged into the 3 h hold."]


# Attach procedure source facts and source-context records.
for r in reactions:
    if r["record_type"] == "specific_experiment" and r["procedure"]["original_text"]["normalized_value"] is None:
        r["procedure"]["original_text"] = F(
            {"location": r["location"]["raw_value"]},
            "See cited table/scheme/procedure evidence; structured steps preserve the reported order.",
            tuple(sorted({int(e.split(":page")[1].split(":")[0]) for e in r["evidence_ids"]})))

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:article_metadata", "document_id": DOC, "kind": "article_metadata",
     "title": TITLE, "doi": "10.1039/c9gc04184b", "evidence_ids": ev(1)},
    {"context_id": f"{DOC}:figure1_xray", "document_id": DOC, "kind": "structure_characterization",
     "description": "Figure 1 X-ray structure of 2a; text states the structure was determined by X-ray crystallography (CCDC 1958998).", "evidence_ids": ev(2)},
    {"context_id": f"{DOC}:characterization_pages4_6", "document_id": DOC, "kind": "characterization",
     "description": "Pages 4-6 contain 1H/13C NMR, HRMS and melting-point characterization attached to the corresponding preparations; these peak lists are not separate synthesis experiments.", "evidence_ids": ev(4, 5, 6)},
    {"context_id": f"{DOC}:si_notice", "document_id": DOC, "kind": "supplementary_information_notice",
     "description": "Page 1 states ESI is available, but no SI file was explicitly supplied or associated in source_mapping.json; no SI facts were imported.", "evidence_ids": ev(1)},
]

st.write(RUN / "submitted.json", bundle)

inventory = {
    "document_id": DOC,
    "title": TITLE,
    "requested_pages": [1, 2, 3, 4, 5, 6],
    "visual_review": "All six registered page images opened at original detail.",
    "record_count": len(reactions),
    "molecule_count": len(molecules),
    "records_by_page": {str(p): [r["reaction_id"] for r in reactions if f"{DOC}:page{p}:text" in r["evidence_ids"]] for p in range(1, 7)},
    "notes": [
        "Table 1 entry 1, Table 2 product 2a and the page-6 specific preparation are one experiment with merged evidence.",
        "Table 2 decomposition/no-reaction cases do not receive invented product participants.",
        "NMR/HRMS/melting-point blocks are source_context, not reaction records.",
        "Pages outside 1-6 were not claimed as visually covered or used to fill missing experimental quantities.",
    ],
}
st.write(RUN / "raw_model/pages_0001_0006_extraction.json", inventory)

stereo_checks = {}
for mid, m in molecules.items():
    smi = st.value(m["smiles"])
    if smi:
        stereo_checks[mid] = st.stereo_inspect(smi)
st.write(RUN / "raw_model/stereo_checks_pages_0001_0006.json", stereo_checks)

print(json.dumps({"molecules": len(molecules), "reactions": len(reactions), "source_context": len(bundle["source_context"])}, ensure_ascii=False))
