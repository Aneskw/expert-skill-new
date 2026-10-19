#!/usr/bin/env python3
"""Build the source-derived v3.1 bundle for PDF pages 1--6.

This is a reproducibility artifact.  It records machine interpretations only;
standard merge will quarantine every item as needs_review.
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
DATASET = json.loads((RUN / "dataset.json").read_text())
DOC = "doc_b27cb890048f595e"
E = {p: f"{DOC}:page{p}:text" for p in range(1, 10)}
I = {p: f"{DOC}:page{p}:image" for p in range(1, 7)}


def f(raw=None, normalized=None, ev=(), method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    if normalized is None and raw is not None:
        normalized = raw
    return st.fact(raw, normalized, ev, method=method, derivation=derivation,
                   status=status, verification=verification, alternatives=alternatives)


def inferred(normalized, ev, derivation):
    return f(None, normalized, ev, method="inference", derivation=derivation)


def derived(normalized, ev, derivation):
    return f(None, normalized, ev, method="structure_derivation", derivation=derivation)


bundle = st.empty(DATASET["documents"])
molecules = {}


def molecule(mid, label, name=None, ev=(E[1],), smiles=None,
             chirality="unresolved", geometry="unresolved", representation=None,
             stereo_note=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    display = name or label
    m["names"]["original"] = f(display, display, ev)
    m["names"]["english"] = f(display, display, ev)
    m["names"]["abbreviations"] = f()
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    m["identity_note"] = f(stereo_note, stereo_note, ev) if stereo_note else f()
    m["external_properties"] = []
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = f(smiles, smiles, ev)
        atom_centers = diag["atom_centers"] or None
        m["stereochemistry"]["atom_centers"] = derived(
            atom_centers, ev,
            "RDKit graph diagnostic on the evidence-linked machine candidate; not source verification."
        ) if atom_centers else f()
        m["stereochemistry"]["double_bonds"] = derived(
            diag["double_bonds"], ev,
            "RDKit graph diagnostic; STEREONONE is not interpreted as an author assignment."
        ) if diag["double_bonds"] else f()
        m["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], ev,
            "Canonical isomeric SMILES calculated from the machine structure candidate."
        )
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], ev,
            "SHA-256 key of the canonical isomeric graph; not a sample identity."
        )
        if representation is None:
            representation = "partial" if any(x["specified"] == "Unspecified" for x in diag["potential_stereo"]) else "encoded"
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f()
        m["stereochemistry"]["atom_centers"] = f()
        m["stereochemistry"]["double_bonds"] = f()
        m["stereochemistry"]["isomeric_smiles"] = f()
        m["stereochemistry"]["identity_key"] = f()
        representation = representation or "not_encoded"
    m["stereochemistry"]["axial"] = f()
    m["stereochemistry"]["mixture"] = f()
    m["stereochemistry"]["assignment"] = f()
    m["stereochemistry"]["chirality_status"] = inferred(
        chirality, ev, "Status assigned from the inspected source drawing/name; absence of stereo marks is not treated as achiral."
    )
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, ev, "E/Z applicability/status assigned from the inspected source drawing/name."
    )
    m["stereochemistry"]["representation_status"] = inferred(
        representation, ev, "Comparison of source stereo content with the retained machine representation."
    )
    molecules[mid] = m
    return mid


# Shared exact small molecules.  Graph diagnostics are generated, but the
# candidates remain unverified after machine merge.
simple = {
    "csf": ("CsF", "cesium fluoride", "[Cs+].[F-]"),
    "acetonitrile": ("CH3CN", "acetonitrile", "CC#N"),
    "dcm": ("CH2Cl2", "dichloromethane", "ClCCl"),
    "thf": ("THF", "tetrahydrofuran", "C1CCOC1"),
    "chloroform": ("CHCl3", "chloroform", "ClC(Cl)Cl"),
    "ethyl_acetate": ("EtOAc", "ethyl acetate", "CCOC(C)=O"),
    "methanol": ("MeOH", "methanol", "CO"),
    "dmf": ("DMF", "N,N-dimethylformamide", "CN(C)C=O"),
    "dmso": ("DMSO", "dimethyl sulfoxide", "CS(C)=O"),
    "tbuok": ("t-BuOK", "potassium tert-butoxide", "CC(C)(C)[O-].[K+]"),
    "naoh": ("NaOH", "sodium hydroxide", "[Na+].[OH-]"),
    "koh": ("KOH", "potassium hydroxide", "[K+].[OH-]"),
    "k2co3": ("K2CO3", "potassium carbonate", "[K+].[K+].[O-]C([O-])=O"),
    "cs2co3": ("Cs2CO3", "cesium carbonate", "[Cs+].[Cs+].[O-]C([O-])=O"),
    "dmap": ("DMAP", "4-(dimethylamino)pyridine", "CN(C)c1ccncc1"),
    "dabco": ("DABCO", "1,4-diazabicyclo[2.2.2]octane", "C1CN2CCN1CC2"),
    "indole": ("indole", "1H-indole", "c1ccc2[nH]ccc2c1"),
    "n_methylindole": ("1-methylindole", "1-methyl-1H-indole", "Cn1ccc2ccccc21"),
    "three_methylindole": ("3-methylindole", "3-methyl-1H-indole", "Cc1c[nH]c2ccccc12"),
    "two_methoxythiophene": ("2-methoxythiophene", "2-methoxythiophene", "COc1cccs1"),
    "n_methylpyrrole": ("1-methylpyrrole", "1-methyl-1H-pyrrole", "Cn1cccc1"),
    "triethylamine": ("Et3N", "triethylamine", "CCN(CC)CC"),
    "boc2o": ("(Boc)2O", "di-tert-butyl dicarbonate", "CC(C)(C)OC(=O)OC(=O)OC(C)(C)C"),
    "nabh3cn": ("NaBH3CN", "sodium cyanoborohydride", "[Na+].[BH3-]C#N"),
    "hcl": ("HCl", "hydrogen chloride", "Cl"),
    "water": ("water", "water", "O"),
    "brine": ("saturated brine", "saturated brine", None),
    "sodium_sulfate": ("Na2SO4", "sodium sulfate", "[Na+].[Na+].[O-]S(=O)(=O)[O-]"),
    "silica": ("silica gel", "silica gel", None),
    "acetone": ("acetone", "acetone", "CC(C)=O"),
    "hexane": ("hexane", "hexane", "CCCCCC"),
}
simple_ev = {
    "csf":(E[3],I[3]), "acetonitrile":(E[1],I[1],E[3],I[3],E[5],I[5]),
    "dcm":(E[1],I[1],E[3],I[3],E[4],I[4],E[6],I[6]), "thf":(E[3],I[3],E[4],I[4]),
    "chloroform":(E[3],I[3]), "ethyl_acetate":(E[3],I[3],E[5],I[5],E[6],I[6]),
    "methanol":(E[3],I[3],E[4],I[4]), "dmf":(E[3],I[3]), "dmso":(E[2],I[2],E[3],I[3]),
    "tbuok":(E[3],I[3]), "naoh":(E[3],I[3]), "koh":(E[3],I[3]), "k2co3":(E[3],I[3]),
    "cs2co3":(E[1],I[1],E[2],I[2],E[3],I[3],E[5],I[5]), "dmap":(E[3],I[3],E[4],I[4]),
    "dabco":(E[3],I[3]), "indole":(E[3],I[3],E[4],I[4],E[6],I[6]),
    "n_methylindole":(E[3],I[3],E[4],I[4]), "three_methylindole":(E[3],I[3],E[4],I[4]),
    "two_methoxythiophene":(E[3],I[3],E[4],I[4]), "n_methylpyrrole":(E[3],I[3],E[4],I[4]),
    "triethylamine":(E[4],I[4]), "boc2o":(E[4],I[4]), "nabh3cn":(E[3],I[3],E[4],I[4]),
    "hcl":(E[3],I[3],E[4],I[4]), "water":(E[2],I[2],E[5],I[5],E[6],I[6]),
    "brine":(E[5],I[5]), "sodium_sulfate":(E[5],I[5],E[6],I[6]),
    "silica":(E[3],I[3],E[5],I[5],E[6],I[6]), "acetone":(E[5],I[5],E[6],I[6]),
    "hexane":(E[3],I[3],E[6],I[6]),
}
for mid, (lab, nm, smi) in simple.items():
    molecule(mid, lab, nm, simple_ev[mid], smi,
             chirality="achiral" if smi else "not_applicable",
             geometry="not_applicable", representation="encoded" if smi else "not_encoded")
# The graph checker flags bridgehead-nitrogen stereo possibilities for DABCO;
# do not overrule that diagnostic with an achiral/fully encoded claim.
molecules["dabco"]["stereochemistry"]["chirality_status"] = inferred(
    "unresolved", (E[3], I[3]), "Toolkit potential-stereo diagnostic is not sufficient to assign or exclude stereochemistry."
)
molecules["dabco"]["stereochemistry"]["representation_status"] = inferred(
    "partial", (E[3], I[3]), "Machine graph retains connectivity but no source stereochemical assignment."
)

# Source-named reagents/complexes whose complete structure is deliberately not guessed.
for mid, lab, nm, ev in [
    ("pph3aucl", "PPh3AuCl", "chlorotriphenylphosphinegold(I)", (E[4], I[4])),
    ("agotf", "AgOTf", "silver triflate", (E[4], I[4])),
    ("dibal_h", "DIBAL-H", "diisobutylaluminium hydride", (E[4], I[4])),
    ("in_otf3", "In(OTf)3", "indium(III) triflate", (E[4], I[4])),
    ("five_bromoindole", "5-bromoindole", "5-bromo-1H-indole", (E[4], I[4])),
    ("six_bromoindole", "6-bromoindole", "6-bromo-1H-indole", (E[4], I[4])),
    ("tetrahydrocyclopentaindole", "tetrahydrocyclopenta[g]indole", "1,6,7,8-tetrahydrocyclopenta[g]indole", (E[4], I[4])),
]:
    molecule(mid, lab, nm, ev, None, chirality="unresolved", geometry="not_applicable",
             stereo_note="Complete source structure was not encoded with sufficient confidence.")

# Scope substrates 1a--t: only 1a gets a machine candidate.  All drawings have
# an unassigned tetrahedral center; the terminal allene is not assigned axial
# configuration and is not silently called racemic.
substrate_labels = {c: f"1{c}" for c in "abcdefghijklmnopqrst"}
molecule("1a", "1a", "ethyl 2-(2-nitro-1-phenylethyl)buta-2,3-dienoate",
         (E[2], E[3], I[3], E[5], I[5]),
         "C=C=C(C(=O)OCC)C(c1ccccc1)C[N+](=O)[O-]",
         chirality="unspecified", geometry="not_applicable", representation="partial",
         stereo_note="One tetrahedral center is present but no configuration is assigned; the terminal allene drawing carries no axial assignment.")
for c in "bcdefghijklmnopqrst":
    molecule(f"1{c}", f"1{c}", f"alpha-nitroethylallenic ester 1{c}",
             (E[3], I[3]), None, chirality="unspecified", geometry="not_applicable",
             stereo_note="Scope drawing identifies the substituent pattern, but a complete stereochemical graph was not encoded; the tetrahedral center is unassigned.")

product_names = {
    "a":"Ethyl 2-hydroxy-2-methyl-5-oxo-4-phenyl-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "b":"Ethyl 2-hydroxy-4-(4-methoxyphenyl)-2-methyl-5-oxo-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "c":"Ethyl 2-hydroxy-2-methyl-4-(4-(methylthio)phenyl)-5-oxo-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "d":"Ethyl 2-hydroxy-2-methyl-5-oxo-4-(p-tolyl)-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "e":"Ethyl 4-([1,1'-biphenyl]-4-yl)-2-hydroxy-2-methyl-5-oxo-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "f":"Ethyl 4-(4-bromophenyl)-2-hydroxy-2-methyl-5-oxo-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "g":"Ethyl 4-(4-chlorophenyl)-2-hydroxy-2-methyl-5-oxo-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "h":"Ethyl 4-(4-fluorophenyl)-2-hydroxy-2-methyl-5-oxo-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "i":"Ethyl 2-hydroxy-4-(3-methoxyphenyl)-2-methyl-5-oxo-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "j":"Ethyl 4-(3-bromophenyl)-2-hydroxy-2-methyl-5-oxo-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "k":"Ethyl 4-(3-fluorophenyl)-2-hydroxy-2-methyl-5-oxo-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "l":"Ethyl 4-(3,4-dimethoxyphenyl)-2-hydroxy-2-methyl-5-oxo-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "m":"Ethyl 2-hydroxy-2-methyl-5-oxo-4-(3,4,5-trimethoxyphenyl)-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "n":"Ethyl 2-hydroxy-2-methyl-4-(naphthalen-2-yl)-5-oxo-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "o":"Ethyl 2-hydroxy-2-methyl-5-oxo-4-(thiophen-2-yl)-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "p":"Ethyl 4-(furan-2-yl)-2-hydroxy-2-methyl-5-oxo-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "q":"Ethyl 4-(furan-3-yl)-2-hydroxy-2-methyl-5-oxo-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "r":"Methyl 2-hydroxy-2-methyl-5-oxo-4-phenyl-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "s":"Benzyl 2-hydroxy-2-methyl-5-oxo-4-phenyl-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "t":"tert-Butyl 2-hydroxy-2-methyl-5-oxo-4-phenyl-2,5-dihydro-1H-pyrrole-3-carboxylate",
}
aryl = {
    "a":"c2ccccc2", "b":"c2ccc(OC)cc2", "c":"c2ccc(SC)cc2", "d":"c2ccc(C)cc2",
    "e":"c2ccc(-c3ccccc3)cc2", "f":"c2ccc(Br)cc2", "g":"c2ccc(Cl)cc2", "h":"c2ccc(F)cc2",
    "i":"c2cc(OC)ccc2", "j":"c2cc(Br)ccc2", "k":"c2cc(F)ccc2",
    "l":"c2cc(OC)c(OC)cc2", "m":"c2cc(OC)c(OC)c(OC)c2", "n":"c2ccc3ccccc3c2",
    "o":"c2cccs2", "p":"c2ccco2", "q":"c2ccoc2",
}
product_smiles = {c: f"CCOC(=O)C1=C({a})C(=O)NC1(C)O" for c, a in aryl.items()}
product_smiles.update({
    "r":"COC(=O)C1=C(c2ccccc2)C(=O)NC1(C)O",
    "s":"O=C(OCc2ccccc2)C1=C(c2ccccc2)C(=O)NC1(C)O",
    "t":"CC(C)(C)OC(=O)C1=C(c2ccccc2)C(=O)NC1(C)O",
})
for c in "abcdefghijklmnopqrst":
    ev = (E[3], I[3], E[5] if c <= "i" else E[6], I[5] if c <= "i" else I[6])
    molecule(f"2{c}", f"2{c}", product_names[c], ev, product_smiles[c],
             chirality="unspecified", geometry="not_applicable", representation="partial",
             stereo_note="The hydroxyl/methyl-bearing tetrahedral center is stereogenic in the graph, but the paper gives no configuration or mixture assignment.")

# Heteroarylation inputs and products.  Products retain source identity and
# unassigned stereogenic centers without invented configurations.
hetero_inputs = {
    "3a":"n_methylindole", "3b":"indole", "3c":"five_bromoindole", "3d":"six_bromoindole",
    "3e":"tetrahydrocyclopentaindole", "3f":"three_methylindole",
    "3g":"two_methoxythiophene", "3h":"n_methylpyrrole",
}
for label, base_mid in hetero_inputs.items():
    # A separate literature label occurrence is represented by the shared
    # concrete molecule identity, not a duplicate molecule object.
    pass
product4_names = {
    "a":"Ethyl 2-methyl-2-(1-methyl-indol-3-yl)-5-oxo-4-phenyl-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "b":"Ethyl 2-(1H-indol-3-yl)-2-methyl-5-oxo-4-phenyl-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "c":"Ethyl 2-(5-bromo-1H-indol-3-yl)-2-methyl-5-oxo-4-phenyl-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "d":"Ethyl 2-(6-bromo-1H-indol-3-yl)-2-methyl-5-oxo-4-phenyl-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "e":"Ethyl 2-methyl-5-oxo-4-phenyl-2-(1,6,7,8-tetrahydrocyclopenta[g]indol-3-yl)-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "f":"Ethyl 2-methyl-2-(3-methyl-1H-indol-2-yl)-5-oxo-4-phenyl-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "g":"Ethyl 2-(5-methoxythiophen-2-yl)-2-methyl-5-oxo-4-phenyl-2,5-dihydro-1H-pyrrole-3-carboxylate",
    "h":"Ethyl 4-(4-methoxyphenyl)-1',2-dimethyl-5-oxo-2,5-dihydro-1H,1'H-[2,2'-bipyrrole]-3-carboxylate",
}
for c in "abcdefgh":
    molecule(f"4{c}", f"4{c}", product4_names[c], (E[4], I[4], E[6] if c <= "c" else E[7]),
             None, chirality="unspecified", geometry="not_applicable", representation="not_encoded",
             stereo_note="The new tetrahedral center is drawn without absolute/relative assignment; no racemic composition is inferred.")

for mid, label, name, chi, geom, mix in [
    ("5", "5", "bis-Boc-protected derivative 5", "unspecified", "not_applicable", None),
    ("6", "6", "hemiaminal 6", "mixture", "not_applicable", "diastereomeric mixture, 5:1"),
    ("7", "7", "bisindole derivative 7", "mixture", "not_applicable", "diastereomeric mixture, 3:1"),
    ("8", "8", "Ethyl 2-methyl-5-oxo-4-phenyl-2,5-dihydro-1H-pyrrole-3-carboxylate", "unspecified", "not_applicable", None),
    ("9", "9", "Ethyl 2-methoxy-2-methyl-5-oxo-4-phenyl-2,5-dihydro-1H-pyrrole-3-carboxylate", "unspecified", "not_applicable", None),
    ("intermediate_f", "Intermediate F", "(Z)-Ethyl 5-(hydroxyimino)-2-methyl-4-phenyl-4,5-dihydrofuran-3-carboxylate", "partially_specified", "specified", None),
]:
    molecule(mid, label, name, (E[4], I[4], E[5], I[5], E[7], E[8]), None,
             chirality=chi, geometry=geom, representation="not_encoded",
             stereo_note=(mix or "Source stereochemical content is retained textually; a complete isomeric graph was not encoded."))
    if mix:
        molecules[mid]["stereochemistry"]["mixture"] = f(mix, {"type":"diastereomeric_mixture", "reported_ratio":mix.rsplit(" ",1)[-1]}, (E[4], I[4]))
    if mid == "intermediate_f":
        molecules[mid]["stereochemistry"]["double_bonds"] = f("(Z)-hydroxyimino", [{"element":"C=N-OH", "configuration":"Z"}], (E[5], I[5], E[8]))
        molecules[mid]["stereochemistry"]["assignment"] = f("(Z) in the reported name", {"type":"geometric", "configuration":"Z"}, (E[8],))

# Generic molecules for the present-work schemes and six prior-work arrows.
for mid, label in [
    ("generic_1", "alpha-nitroethylallenic ester 1"), ("generic_2", "5-hydroxy-3-pyrrolin-2-one 2"),
    ("generic_3", "electron-rich heteroarene 3"), ("generic_4", "heteroaryl-lactam 4"),
]:
    molecule(mid, label, label, (E[1], I[1], E[2], I[2], E[4], I[4]), None,
             chirality="unresolved", geometry="unresolved", representation="not_encoded",
             stereo_note="Generic R/Ar-variable drawing; retained without a false concrete SMILES.")
for n in range(1, 7):
    molecule(f"lit_r{n}", f"Scheme 1 literature reactant {n}", ev=(E[2], I[2]),
             chirality="unresolved", geometry="unresolved", representation="not_encoded",
             stereo_note="Generic prior-work structure with R variables.")
    molecule(f"lit_p{n}", f"Scheme 1 literature product {n}", ev=(E[2], I[2]),
             chirality="unresolved", geometry="unresolved", representation="not_encoded",
             stereo_note="Generic prior-work hydroxypyrrolinone structure with R variables.")
for mid, label in [
    ("tbuocl", "t-BuOCl"), ("iodine", "I2"), ("cubr2", "CuBr2"), ("can", "CAN"),
    ("fecl3", "FeCl3"), ("acoh", "AcOH"), ("amine_generic", "R4NH2"),
]:
    molecule(mid, label, label, (E[2], I[2]), None, chirality="not_applicable", geometry="not_applicable")


def raw_quantity(value):
    if not isinstance(value, dict):
        return str(value)
    if "number" in value:
        prefix = {"lt":"<", "le":"≤", "gt":">", "ge":"≥", "approx":"~"}.get(value.get("qualifier"), "")
        return f"{prefix}{value['number']} {value.get('unit','')}".strip()
    if "minimum" in value or "maximum" in value:
        return f"{value.get('minimum','?')}–{value.get('maximum','?')} {value.get('unit','')}".strip()
    return json.dumps(value, ensure_ascii=False)


def charge_fact(value, ev):
    return f(raw_quantity(value), value, ev)


def participant(pid, mid, role, ev, *, stage="reaction", mass=None, amount=None,
                volume=None, equiv=None, mol_percent=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage)
    p["product_role"] = f(product_role, product_role, ev) if product_role else f()
    values = {"mass":mass, "amount":amount, "volume":volume,
              "equivalents":equiv, "mol_percent":mol_percent}
    for key, value in values.items():
        p["charges"][key] = charge_fact(value, ev) if value is not None else f()
    p["equivalence_reference"] = f(reference, reference, ev) if reference else f()
    return p


def step(sid, ev, pids, operation, temperature=None, time=None, stage="reaction",
         atmosphere=None, pressure=None):
    s = st.new_step(sid, ev, stage)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, ev)
    s["temperature"] = f(raw_quantity(temperature), temperature, ev) if temperature is not None else f()
    s["time"] = f(raw_quantity(time), time, ev) if time is not None else f()
    s["atmosphere"] = f(atmosphere, atmosphere, ev) if atmosphere else f()
    s["pressure"] = f(pressure, pressure, ev) if pressure else f()
    return s


def bond_class(ev, formed=(), h_gains=(), limitations=()):
    changes = []
    for a, b in formed:
        changes.append({"change_id":f"b{len(changes)+1}", "atoms":[{"atom_map":None,"element":a},{"atom_map":None,"element":b}],
                        "before_order":0, "after_order":1, "change_type":"formed"})
    for element in h_gains:
        changes.append({"change_id":f"b{len(changes)+1}", "atoms":[{"atom_map":None,"element":element},{"atom_map":None,"element":"H"}],
                        "before_order":0, "after_order":1, "change_type":"hydrogen_gain", "hydrogen_count":1})
    labels = {}
    for ch in changes:
        lab = st.bond_label(ch["atoms"][0]["element"], ch["atoms"][1]["element"], ch["after_order"])
        labels.setdefault(lab["code"], {**lab, "supporting_change_ids":[]})["supporting_change_ids"].append(ch["change_id"])
    deriv = "Partial net-bond interpretation from the inspected reactant/product drawings; no atom mapping is asserted."
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": f(),
        "bond_changes": [derived(ch, ev, deriv) for ch in changes],
        "labels": derived(list(labels.values()), ev, deriv) if changes else inferred([], ev, "No formed-bond label asserted."),
        "coverage": inferred("partial" if changes else "unresolved", ev, deriv if changes else "Source structures are too generic/ambiguous for safe classification."),
        "limitations": list(limitations) or ["No atom-mapped reaction SMILES; leaving-group, atom-source, and rearrangement changes may be incomplete."],
        "model_votes": [],
    }


def reaction(rid, ev, record_type, categories, location, name, transformation,
             participants, steps=(), outcomes=(), classification=None, notes=(), unresolved=(),
             general_scheme_ids=(), solvent=None, volume=None, scale=None, procedure=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, ev)
    r["location"] = f(location, location, ev)
    r["general_scheme_ids"] = list(general_scheme_ids)
    r["reaction_name"] = f(name, name, ev)
    r["transformation_type"] = inferred(transformation, ev, "Machine interpretation of the source-described transformation.")
    r["participants"] = list(participants)
    r["reaction_scale"] = f(f"{raw_quantity(scale['amount'])} of {scale['reference_molecule_id']}", scale, ev) if scale else f()
    r["solvent_conditions"]["composition"] = f(solvent, solvent, ev) if solvent else f()
    r["solvent_conditions"]["mixture_ratio"] = f()
    r["solvent_conditions"]["total_volume"] = f(raw_quantity(volume), volume, ev) if volume else f()
    r["solvent_conditions"]["concentration"] = f()
    r["solvent_conditions"]["concentration_basis"] = f()
    r["steps"] = list(steps)
    r["procedure"]["original_text"] = f(procedure, procedure, ev) if procedure else f()
    r["procedure"]["general_procedure_reference"] = f()
    r["procedure"]["modifications"] = f()
    r["procedure"]["addition_order"] = f()
    r["procedure"]["pretreatment"] = f()
    r["outcomes"] = list(outcomes)
    r["classification"] = classification or bond_class(ev)
    r["notes"] = [f(n, n, ev) for n in notes]
    r["review"] = {"status":"needs_review", "method":"model", "reviewer":None,
                   "checked_fields":[], "unresolved":list(unresolved)}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r["reaction_id"]


def outcome(oid, kind, target, raw, normalized, ev, basis=None, yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, [target] if target else [], ev)
    o["value"] = f(raw, normalized, ev)
    o["product_role"] = f("main product", "main_product", ev) if kind in ("yield", "product_amount") else f()
    o["yield_type"] = f(yield_type, yield_type, ev) if yield_type else f()
    o["measurement_basis"] = f(basis, basis, ev) if basis else f()
    o["scope"] = f(scope, scope, ev) if scope else f()
    return o


record_ids_by_page = {p: [] for p in range(1, 7)}
def link(rid, *pages):
    for p in pages:
        record_ids_by_page[p].append(rid)


# Six prior-work arrows from Scheme 1 (two each in a--c).
literature = [
    (1, ["tbuocl", "iodine"], "MeCN", "halolactamization of a homoallylic amide"),
    (2, ["cubr2"], "THF/H2O (1:1)", "Cu-mediated haloacetamization of an allenamide"),
    (3, ["can"], "MeCN", "oxidative cyclization of an enamide"),
    (4, ["fecl3"], "CH2Cl2", "FeCl3-catalyzed intramolecular enaminic addition"),
    (5, [], "DMSO/H2O; 80 °C", "intermolecular nitroolefin/1,3-diketone tandem reaction"),
    (6, ["amine_generic", "acoh"], "toluene; reflux", "amine addition/cyclization of an MBH adduct"),
]
for n, reagents, solvent_desc, desc in literature:
    ev = (E[2], I[2])
    ps = [participant("r", f"lit_r{n}", "reactant", ev), participant("p", f"lit_p{n}", "product", ev, product_role="main product")]
    ps += [participant(f"x{i}", m, "reagent", ev) for i, m in enumerate(reagents, 1)]
    rid = reaction(f"scheme1_literature_arrow_{n}", ev, "literature_example", [],
                   {"page":2, "scheme":"Scheme 1", "arrow":n}, desc, desc, ps,
                   [step("reaction", ev, [p["participant_id"] for p in ps], solvent_desc)],
                   classification=bond_class(ev),
                   unresolved=("Generic R-group prior-work structures and yields are not reported in this review scheme.",))
    link(rid, 2)

# Present-work generic schemes, also shown in the page-1 graphical abstract.
ev12 = (E[1], I[1], E[2], I[2], E[3], I[3], E[5], I[5])
ps = [participant("s", "generic_1", "reactant", ev12), participant("p", "generic_2", "product", ev12, product_role="main product"),
      participant("base", "cs2co3", "reagent", ev12, equiv={"number":1,"unit":"equiv"}, reference="generic_1"),
      participant("solv", "acetonitrile", "solvent", ev12)]
generic_cyclization = reaction("general_cyclization_1_to_2", ev12, "generic_scheme", [],
    {"pages":[1,2,3,5], "schemes":["graphical abstract","Scheme 1d","Scheme 2","general procedure"]},
    "General procedure for the synthesis of 5-hydroxy-3-pyrrolin-2-ones 2a-t",
    "base-mediated intramolecular cyclization/rearrangement", ps,
    [step("reaction", ev12, ["s","base","solv","p"], "stir", {"number":25,"unit":"°C"}, {"number":2,"unit":"h"})],
    classification=bond_class(ev12, formed=(("C","N"),), h_gains=("N","O"),
        limitations=("Partial net interpretation only; the rearrangement/hydrolysis atom mapping and oxygen/hydrogen sources are unresolved.",)),
    notes=("The graphical abstract summarizes 20 examples in 61-88% yield; this is a scope summary, not a separate experiment.",),
    unresolved=("Generic Ar/R structures are variable; full atom mapping is not established.",),
    solvent="CH3CN", volume={"number":2,"unit":"mL"})
link(generic_cyclization, 1, 2, 3, 5)

ev146 = (E[1], I[1], E[2], I[2], E[4], I[4], E[6], I[6])
ps = [participant("lactam", "generic_2", "reactant", ev146), participant("het", "generic_3", "reactant", ev146),
      participant("p", "generic_4", "product", ev146, product_role="main product"),
      participant("au", "pph3aucl", "catalyst", ev146, mol_percent={"number":10,"unit":"mol%"}, reference="generic_2"),
      participant("ag", "agotf", "catalyst", ev146, mol_percent={"number":10,"unit":"mol%"}, reference="generic_2"),
      participant("solv", "dcm", "solvent", ev146)]
generic_hetero = reaction("general_heteroarylation_2_plus_3_to_4", ev146, "generic_scheme", [],
    {"pages":[1,2,4,6], "schemes":["graphical abstract","Scheme 1d","Scheme 3","general procedure"]},
    "General procedure for the heteroarylation of 2-hydroxy-2,5-dihydro-1H-pyrrole-3-carboxylates 4a-h",
    "Au/Ag-catalyzed heteroarylation", ps,
    [step("reaction", ev146, [p["participant_id"] for p in ps], "stir", {"number":25,"unit":"°C"}, {"number":30,"unit":"min"}, atmosphere="air")],
    classification=bond_class(ev146, formed=(("C","C"),),
        limitations=("New heteroaryl-to-lactam C-C connection only; C-H/O-H losses and exact atom mapping are unresolved.",)),
    notes=("The graphical abstract summarizes 8 examples in 66-86% yield; this is a scope summary, not a separate experiment.",),
    unresolved=("Generic heteroarene and lactam structures are variable.",), solvent="CH2Cl2", volume={"number":2,"unit":"mL"})
link(generic_hetero, 1, 2, 4, 6)

# Table 1: all 17 rows are independent experiments.
table = [
    (1,"csf","acetonitrile",12,66,None,25),(2,"csf","dcm",12,15,None,25),(3,"csf","thf",12,"<5",None,25),
    (4,"csf","chloroform",12,"<5",None,25),(5,"csf","ethyl_acetate",12,11,None,25),(6,"csf","methanol",12,"<5",None,25),
    (7,"csf","dmf",12,16,None,25),(8,"csf","dmso",12,12,None,25),(9,"tbuok","acetonitrile",12,15,None,25),
    (10,"naoh","acetonitrile",2,40,None,25),(11,"koh","acetonitrile",2,65,None,25),(12,"k2co3","acetonitrile",2,61,None,25),
    (13,"cs2co3","acetonitrile",2,86,None,25),(14,"dmap","acetonitrile",12,55,None,25),(15,"dabco","acetonitrile",12,75,None,25),
    (16,"cs2co3","acetonitrile",12,66,0.10,25),(17,"cs2co3","acetonitrile",2,85,None,70),
]
for entry, base, solv, hours, y, special_amount, temp in table:
    ev = (E[3], I[3])
    base_amt = special_amount if special_amount is not None else 0.20
    ps = [participant("s", "1a", "reactant", ev, amount={"number":0.20,"unit":"mmol"}),
          participant("p", "2a", "product", ev, product_role="main product"),
          participant("base", base, "reagent", ev, amount={"number":base_amt,"unit":"mmol"},
                      equiv={"number":base_amt/0.20,"unit":"equiv"}, reference="1a"),
          participant("solv", solv, "solvent", ev, volume={"number":2.0,"unit":"mL"})]
    ynorm = {"number":y,"unit":"%"} if isinstance(y,int) else {"number":5,"unit":"%","qualifier":"lt","raw_text":"<5%"}
    outs = [outcome("yield_2a", "yield", "2a", f"{y}%" if isinstance(y,int) else "<5%", ynorm, ev,
                    basis="isolated after silica gel column chromatography", yield_type="isolated")]
    rid = reaction(f"table1_entry_{entry:02d}", ev, "specific_experiment", ["optimization"],
        {"page":3,"table":"Table 1","entry":entry}, "Optimization of the reaction conditions",
        "base-mediated intramolecular cyclization/rearrangement", ps,
        [step("reaction", ev, ["s","base","solv","p"], "stir", {"number":temp,"unit":"°C"}, {"number":hours,"unit":"h"})],
        outs, classification=bond_class(ev, formed=(("C","N"),), h_gains=("N","O"),
            limitations=("Partial net interpretation; rearrangement/hydrolysis atom mapping is unresolved.",)),
        unresolved=("The table does not report workup details for individual optimization rows.",),
        general_scheme_ids=(generic_cyclization,), solvent=simple[solv][0], volume={"number":2.0,"unit":"mL"},
        scale={"reference_molecule_id":"1a","amount":{"number":0.20,"unit":"mmol"}})
    link(rid, 3)

# Scheme 2 scope.
scope_yields = dict(zip("abcdefghijklmnopqrst", [86,82,85,74,86,71,69,81,81,74,61,71,82,88,61,65,64,71,68,67]))
for c, y in scope_yields.items():
    ev = (E[3], I[3], E[5] if c <= "i" else E[6], I[5] if c <= "i" else I[6])
    ps = [participant("s", f"1{c}", "reactant", ev, amount={"number":0.20,"unit":"mmol"}),
          participant("p", f"2{c}", "product", ev, product_role="main product"),
          participant("base", "cs2co3", "reagent", ev, amount={"number":0.20,"unit":"mmol"}, equiv={"number":1,"unit":"equiv"}, reference=f"1{c}"),
          participant("solv", "acetonitrile", "solvent", ev, volume={"number":2.0,"unit":"mL"})]
    out = outcome(f"yield_2{c}", "yield", f"2{c}", f"{y}%", {"number":y,"unit":"%"}, ev,
                  basis="isolated after silica gel column chromatography", yield_type="isolated")
    rid = reaction(f"scheme2_scope_2{c}", ev, "specific_experiment", ["substrate_scope"],
        {"page":3,"scheme":"Scheme 2","product":f"2{c}"}, "Synthesis of 5-hydroxy-3-pyrrolin-2-ones",
        "base-mediated intramolecular cyclization/rearrangement", ps,
        [step("reaction", ev, ["s","base","solv","p"], "stir", {"number":25,"unit":"°C"}, {"number":2,"unit":"h"})],
        [out], classification=bond_class(ev, formed=(("C","N"),), h_gains=("N","O"),
            limitations=("Partial net interpretation; rearrangement/hydrolysis atom mapping is unresolved.",)),
        unresolved=("Starting-material preparation is cited to prior work; no associated SI was supplied.",),
        general_scheme_ids=(generic_cyclization,), solvent="CH3CN", volume={"number":2.0,"unit":"mL"},
        scale={"reference_molecule_id":f"1{c}","amount":{"number":0.20,"unit":"mmol"}})
    link(rid, 3, 5 if c <= "i" else 6)

# Distinct gram-scale experiment for 1a -> 2a.
ev = (E[3], I[3], E[5], I[5])
ps = [participant("s", "1a", "reactant", ev, mass={"number":1.0,"unit":"g"}, amount={"number":3.82,"unit":"mmol"}),
      participant("p", "2a", "product", ev, product_role="main product"),
      participant("base", "cs2co3", "reagent", ev, mass={"number":1.247,"unit":"g"}, amount={"number":3.82,"unit":"mmol"}, equiv={"number":1,"unit":"equiv"}, reference="1a"),
      participant("solv", "acetonitrile", "solvent", ev, volume={"number":5,"unit":"mL"}),
      participant("extract", "ethyl_acetate", "workup_agent", ev, stage="workup", volume={"number":20,"unit":"mL"}),
      participant("water", "water", "workup_agent", ev, stage="workup"), participant("brine", "brine", "workup_agent", ev, stage="workup"),
      participant("dry", "sodium_sulfate", "workup_agent", ev, stage="workup"), participant("silica", "silica", "purification_agent", ev, stage="purification"),
      participant("acetone", "acetone", "purification_agent", ev, stage="purification"), participant("dcm_p", "dcm", "purification_agent", ev, stage="purification")]
outs = [outcome("yield_2a", "yield", "2a", "66%", {"number":66,"unit":"%"}, ev, basis="isolated after silica gel column chromatography", yield_type="isolated"),
        outcome("amount_2a", "product_amount", "2a", "660 mg", {"number":660,"unit":"mg"}, ev, basis="isolated product", scope="gram-scale experiment")]
rid = reaction("gram_scale_1a_to_2a", ev, "specific_experiment", ["scale_up"],
    {"pages":[3,5],"scheme":"Scheme 2 footnote b and gram-scale procedure"}, "Gram-scale synthesis of 2a",
    "base-mediated intramolecular cyclization/rearrangement", ps,
    [step("reaction", ev, ["s","base","solv","p"], "stir", {"number":25,"unit":"°C"}, {"number":2,"unit":"h"}),
     step("workup", ev, ["extract","water","brine","dry"], "extract with ethyl acetate (3 x 20 mL) and water; wash with saturated brine; dry and evaporate", stage="workup"),
     step("purification", ev, ["silica","acetone","dcm_p"], "silica-gel column chromatography using acetone/dichloromethane", stage="purification")],
    outs, classification=bond_class(ev, formed=(("C","N"),), h_gains=("N","O"),
        limitations=("Partial net interpretation; rearrangement/hydrolysis atom mapping is unresolved.",)),
    general_scheme_ids=(generic_cyclization,), solvent="CH3CN", volume={"number":5,"unit":"mL"},
    scale={"reference_molecule_id":"1a","amount":{"number":3.82,"unit":"mmol"}})
link(rid, 3, 5)

# Scheme 3 scope.
hetero_yields = dict(zip("abcdefgh", [86,77,67,66,66,80,72,67]))
hetero_mid = {"a":"n_methylindole","b":"indole","c":"five_bromoindole","d":"six_bromoindole","e":"tetrahydrocyclopentaindole","f":"three_methylindole","g":"two_methoxythiophene","h":"n_methylpyrrole"}
for c, y in hetero_yields.items():
    ev = (E[4], I[4], E[6], I[6])
    start = "2b" if c == "h" else "2a"
    ps = [participant("lactam", start, "reactant", ev, amount={"number":0.20,"unit":"mmol"}),
          participant("het", hetero_mid[c], "reactant", ev, amount={"number":0.26,"unit":"mmol"}, equiv={"number":1.3,"unit":"equiv"}, reference=start),
          participant("p", f"4{c}", "product", ev, product_role="main product"),
          participant("au", "pph3aucl", "catalyst", ev, amount={"number":0.02,"unit":"mmol"}, mol_percent={"number":10,"unit":"mol%"}, reference=start),
          participant("ag", "agotf", "catalyst", ev, amount={"number":0.02,"unit":"mmol"}, mol_percent={"number":10,"unit":"mol%"}, reference=start),
          participant("solv", "dcm", "solvent", ev, volume={"number":2.0,"unit":"mL"}),
          participant("extract", "ethyl_acetate", "workup_agent", ev, stage="workup"), participant("water", "water", "workup_agent", ev, stage="workup"),
          participant("dry", "sodium_sulfate", "workup_agent", ev, stage="workup"), participant("silica", "silica", "purification_agent", ev, stage="purification"),
          participant("hexane", "hexane", "purification_agent", ev, stage="purification")]
    outs = [outcome(f"yield_4{c}", "yield", f"4{c}", f"{y}%", {"number":y,"unit":"%"}, ev,
                    basis="isolated after silica gel column chromatography", yield_type="isolated")]
    masses = {"a":64,"b":56,"c":59}
    if c in masses:
        outs.append(outcome(f"amount_4{c}", "product_amount", f"4{c}", f"{masses[c]} mg", {"number":masses[c],"unit":"mg"}, ev,
                            basis="isolated product", scope="product characterization block"))
    rid = reaction(f"scheme3_scope_4{c}", ev, "specific_experiment", ["substrate_scope"],
        {"page":4,"scheme":"Scheme 3","product":f"4{c}"}, "Heteroarylation of 5-hydroxy-3-pyrrolin-2-ones",
        "Au/Ag-catalyzed heteroarylation", ps,
        [step("reaction", ev, ["lactam","het","au","ag","solv","p"], "stir", {"number":25,"unit":"°C"}, {"number":30,"unit":"min"}, atmosphere="air"),
         step("workup", ev, ["extract","water","dry"], "evaporate solvent; extract with ethyl acetate and water; dry and evaporate", stage="workup"),
         step("purification", ev, ["silica","extract","hexane"], "silica-gel chromatography, ethyl acetate/hexane 40/60", stage="purification")],
        outs, classification=bond_class(ev, formed=(("C","C"),),
            limitations=("Only the newly drawn heteroaryl-lactam C-C connection is classified; exact atom mapping and proton loss are unresolved.",)),
        unresolved=("Absolute configuration at the new tetrahedral center is not reported.",),
        general_scheme_ids=(generic_hetero,), solvent="CH2Cl2", volume={"number":2.0,"unit":"mL"},
        scale={"reference_molecule_id":start,"amount":{"number":0.20,"unit":"mmol"}})
    link(rid, 4, 6)

# Scheme 4 derivatizations.  Page 7/8 extracted text is used only as context;
# semantic page coverage remains strictly 1--6.
def deriv_record(rid_suffix, source_mid, product_mid, extras, reaction_steps, y, mass, dr=None, class_args=None, name="Synthetic modification"):
    ev = (E[4], I[4], E[7], E[8])
    ps = [participant("s", source_mid, "reactant", ev), participant("p", product_mid, "product", ev, product_role="main product")]
    ps += extras
    outs = [outcome(f"yield_{product_mid}", "yield", product_mid, f"{y}%", {"number":y,"unit":"%"}, ev,
                    basis="isolated after column chromatography", yield_type="isolated")]
    if mass:
        outs.append(outcome(f"amount_{product_mid}", "product_amount", product_mid, f"{mass} mg", {"number":mass,"unit":"mg"}, ev, basis="isolated product"))
    if dr:
        outs.append(outcome(f"dr_{product_mid}", "selectivity", product_mid, dr, {"metric":"dr","ratio":dr,"order":"major/minor unassigned"}, ev,
                            basis=None, scope="diastereomer ratio; configurations and measurement method unreported"))
    rr = reaction(rid_suffix, ev, "specific_experiment", ["derivatization"], {"page":4,"scheme":"Scheme 4"}, name,
        name, ps, reaction_steps, outs,
        classification=bond_class(ev, **(class_args or {})),
        unresolved=(("Diastereomer-ratio measurement method and configurations are unreported." if dr else "Absolute configuration is unreported."),))
    link(rr, 4)
    return rr

ev478 = (E[4], I[4], E[7], E[8])
deriv_record("scheme4_4b_to_5", "4b", "5", [
    participant("boc", "boc2o", "reagent", ev478, mass={"number":218,"unit":"mg"}, amount={"number":1.0,"unit":"mmol"}, equiv={"number":2.5,"unit":"equiv"}, reference="4b"),
    participant("et3n", "triethylamine", "reagent", ev478, mass={"number":40,"unit":"mg"}, amount={"number":0.4,"unit":"mmol"}, equiv={"number":1,"unit":"equiv"}, reference="4b"),
    participant("dmap", "dmap", "catalyst", ev478, mass={"number":24,"unit":"mg"}, amount={"number":0.2,"unit":"mmol"}, equiv={"number":0.5,"unit":"equiv"}, reference="4b"),
    participant("solv", "thf", "solvent", ev478, volume={"number":3,"unit":"mL"})],
    [step("addition", ev478, ["s","boc","et3n","dmap","solv"], "add (Boc)2O slowly", {"number":0,"unit":"°C"}),
     step("hold", ev478, ["s","boc","et3n","dmap","solv","p"], "warm to 25 °C and stir", {"number":25,"unit":"°C"}, {"number":8,"unit":"h"})],
    55, 122, class_args={"formed":(("C","N"),("C","N")), "limitations":("Two N-Boc connections inferred from the doubly protected product; full atom mapping is absent.",)},
    name="N-Boc protection of 4b")

deriv_record("scheme4_5_to_6", "5", "6", [
    participant("red", "dibal_h", "reagent", ev478, volume={"number":0.8,"unit":"mL"}, equiv={"number":4,"unit":"equiv"}, reference="5"),
    participant("solv", "thf", "solvent", ev478, volume={"number":3,"unit":"mL"})],
    [step("reaction", ev478, ["s","red","solv","p"], "add 1.0 M DIBAL-H in hexane and stir", {"number":-78,"unit":"°C"}, {"number":30,"unit":"min"}, atmosphere="N2")],
    73, 82, dr="5:1", class_args={"h_gains":("C",), "limitations":("Net carbonyl reduction represented by C-H gain; C-O bond-order decrease is not atom-mapped.",)},
    name="DIBAL-H reduction of 5")

deriv_record("scheme4_6_to_7", "6", "7", [
    participant("ind", "n_methylindole", "reactant", ev478, mass={"number":39,"unit":"mg"}, amount={"number":0.30,"unit":"mmol"}, equiv={"number":1.5,"unit":"equiv"}, reference="6"),
    participant("cat", "in_otf3", "catalyst", ev478, mass={"number":11,"unit":"mg"}, amount={"number":0.02,"unit":"mmol"}, mol_percent={"number":10,"unit":"mol%"}, reference="6"),
    participant("solv", "dcm", "solvent", ev478, volume={"number":2,"unit":"mL"})],
    [step("reaction", ev478, ["s","ind","cat","solv","p"], "stir", {"number":25,"unit":"°C"}, {"number":15,"unit":"min"}, atmosphere="air")],
    67, 90, dr="3:1", class_args={"formed":(("C","C"),), "limitations":("New indolyl-hemiaminal carbon connection only; water loss/proton changes are unresolved.",)},
    name="In(OTf)3-catalyzed indolylation of 6")

ev48 = (E[4], I[4], E[7], E[8])
deriv_record("scheme4_2a_to_8", "2a", "8", [
    participant("red", "nabh3cn", "reagent", ev48, mass={"number":38,"unit":"mg"}, amount={"number":0.60,"unit":"mmol"}, equiv={"number":3,"unit":"equiv"}, reference="2a"),
    participant("acid", "hcl", "reagent", ev48, equiv={"number":1.5,"unit":"equiv"}, reference="2a"),
    participant("solv", "thf", "solvent", ev48, volume={"number":2,"unit":"mL"})],
    [step("pre_reduction", ev48, ["s","red","solv"], "add NaBH3CN; after 5 min acidify", {"number":25,"unit":"°C"}, {"number":5,"unit":"min"}),
     step("acidified_hold", ev48, ["s","red","acid","solv","p"], "stir after addition of 35% HCl", {"number":25,"unit":"°C"}, {"number":5,"unit":"h"})],
    56, 27, class_args={"h_gains":("C",), "limitations":("Net replacement of C-OH by C-H; the C-O cleavage is not atom-mapped.",)},
    name="Reductive removal of the hydroxyl group of 2a")

deriv_record("scheme4_2a_to_9", "2a", "9", [
    participant("acid", "hcl", "reagent", ev48, equiv={"number":1,"unit":"equiv"}, reference="2a"),
    participant("solv", "methanol", "solvent", ev48, volume={"number":2,"unit":"mL"})],
    [step("reaction", ev48, ["s","acid","solv","p"], "acidify with 35% HCl and stir", {"number":25,"unit":"°C"}, {"number":1,"unit":"h"})],
    70, 39, class_args={"formed":(("C","O"),), "limitations":("Methoxy C-O connection inferred; identity of the retained/displaced oxygen is not atom-mapped.",)},
    name="Acid-catalyzed methoxylation of 2a")

# Proposed pathway is one nonexperimental record.
ev4 = (E[3], I[3], E[4], I[4])
rid = reaction("scheme5_proposed_mechanism", ev4, "proposed_mechanism", [],
    {"pages":[3,4],"scheme":"Scheme 5"}, "Proposed mechanistic pathway",
    "proposed multistep cyclization/rearrangement pathway through A-I", [
        participant("s", "generic_1", "reactant", ev4), participant("base", "cs2co3", "reagent", ev4), participant("p", "generic_2", "product", ev4)],
    classification=bond_class(ev4),
    notes=("Authors propose deprotonation, nitronate conjugate addition, ring contraction/tautomerization, oxime ring opening, hydrolysis, tautomerization, and final cyclization.",),
    unresolved=("Intermediates A-I are proposed rather than isolated experiments; no atom mapping is asserted.",))
link(rid, 3, 4)

# Scheme 6 is a performed mechanistic control.
ev5 = (E[4], I[4], E[5], I[5])
ps = [participant("s", "intermediate_f", "reactant", ev5), participant("p", "2a", "product", ev5, product_role="main product"),
      participant("base", "cs2co3", "reagent", ev5, equiv={"number":1,"unit":"equiv"}, reference="intermediate_f"), participant("solv", "acetonitrile", "solvent", ev5)]
rid = reaction("scheme6_intermediate_f_to_2a", ev5, "specific_experiment", ["mechanistic_study"],
    {"page":5,"scheme":"Scheme 6"}, "Mechanistic investigation: Intermediate F to 2a",
    "base-mediated conversion of oxime intermediate F to hydroxypyrrolinone 2a", ps,
    [step("reaction", ev5, ["s","base","solv","p"], "stir", {"number":25,"unit":"°C"}, {"number":2,"unit":"h"})],
    [outcome("yield_2a", "yield", "2a", "87%", {"number":87,"unit":"%"}, ev5, basis=None, yield_type=None)],
    classification=bond_class(ev5, formed=(("C","N"),), h_gains=("N",),
        limitations=("Only the final ring-closing C-N connection and net N-H gain are represented; ring opening/hydrolysis and atom sources are unresolved.",)),
    unresolved=("Scheme 6 does not state reaction scale, solvent volume, workup, purification, or yield measurement basis.",),
    general_scheme_ids=(generic_cyclization,), solvent="CH3CN")
link(rid, 5)

# Source context: characterization and nonreaction observations are retained
# without generating synthetic records.
bundle["source_context"] = [
    {"context_id":f"{DOC}:figure1_natural_products", "document_id":DOC, "evidence_ids":[E[2],I[2]], "kind":"nonreaction_figure", "description":"Figure 1 shows bioactive natural products containing the hydroxypyrrolinone motif; it is not a synthesis experiment."},
    {"context_id":f"{DOC}:xray_2f", "document_id":DOC, "evidence_ids":[E[2],E[3],I[3]], "kind":"structure_characterization", "description":"Product 2f was confirmed by X-ray analysis; no separate synthesis record is created."},
    {"context_id":f"{DOC}:xray_4e", "document_id":DOC, "evidence_ids":[E[3],E[4],I[4]], "kind":"structure_characterization", "description":"Product 4e was confirmed by X-ray crystallographic analysis; no separate synthesis record is created."},
    {"context_id":f"{DOC}:pages5_6_characterization", "document_id":DOC, "evidence_ids":[E[5],I[5],E[6],I[6]], "kind":"characterization_blocks", "description":"Melting points, Rf, NMR and HRMS blocks for 2a-t and 4a-c are characterization evidence, not independent reactions."},
    {"context_id":f"{DOC}:gamma_substitution_failure", "document_id":DOC, "evidence_ids":[E[3],I[3]], "kind":"qualitative_scope_limit", "description":"Attempts to introduce gamma substitution in the allenoate failed, but substrate identity and experimental conditions are not specified sufficiently for a separate experiment."},
    {"context_id":f"{DOC}:si_notice_unassociated", "document_id":DOC, "evidence_ids":[E[8]], "kind":"source_availability", "description":"The article states that SI is available online, but no SI file was explicitly provided/associated in this workspace."},
]

bundle["molecules"] = list(molecules.values())
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

inventory = {
    "document_id":DOC,
    "page_range":[1,6],
    "chemical_accuracy_verified":False,
    "molecule_count":len(bundle["molecules"]),
    "reaction_count":len(bundle["reactions"]),
    "counts":{
        "literature_examples":6,
        "generic_present_work_schemes":2,
        "table1_optimization_rows":17,
        "scheme2_scope_rows":20,
        "gram_scale_experiments":1,
        "scheme3_scope_rows":8,
        "scheme4_derivatizations":5,
        "proposed_mechanisms":1,
        "mechanistic_controls":1,
    },
    "record_ids_by_page":record_ids_by_page,
    "notes":[
        "Figure 1 and analytical peak lists were kept as source_context.",
        "The qualitative gamma-substitution failure lacks a substrate identity/conditions and was not fabricated into a row.",
        "Later-page extracted text was used only to complete Scheme 4 quantities; semantic coverage is claimed only for pages 1-6.",
    ],
}
(RUN / "raw_model/pages_0001_0006_inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")

coverage_regions = {
    1:[("title/graphical abstract", "Article identity and two present-work generic transformations; graphical scope summaries are not specific experiments."),
       ("abstract and introduction", "Current-work overview plus cited background; no additional specific experiment." )],
    2:[("Figure 1", "Nonreaction natural-product context."),
       ("Scheme 1a-c", "Six prior-work arrows retained as literature examples."),
       ("Scheme 1d", "Two present-work generic transformations, duplicate of the graphical abstract/general schemes."),
       ("optimization/scope prose", "Discussion cross-references Table 1 and Scheme 2; not duplicated as experiments." )],
    3:[("Table 1", "All 17 optimization rows and footnotes inventoried as independent experiments."),
       ("Scheme 2", "All 20 scope products and footnotes inventoried; gram-scale 1a experiment distinguished."),
       ("discussion", "Scheme 3/4 yields and scope limitations cross-checked; no duplicate records." )],
    4:[("Scheme 3", "Eight heteroarylation scope arrows/outputs inventoried."),
       ("Scheme 4", "Five distinct synthetic modifications inventoried, including dr results for 6 and 7."),
       ("Scheme 5", "Proposed mechanism retained as one nonexperimental mechanism record."),
       ("X-ray 4e", "Structure-characterization context, not a synthesis." )],
    5:[("Scheme 6", "Performed Intermediate F to 2a mechanistic control, 87% yield."),
       ("general experimental information", "Analytical/general laboratory context; no standalone reaction."),
       ("general and gram-scale procedures", "Detailed cyclization and gram-scale charges/workup linked to existing generic/specific records."),
       ("2a-2i characterization", "Product masses/yields and analytical blocks linked to Scheme 2; peak lists are not reactions." )],
    6:[("2i-2t characterization", "Product masses/yields and analytical blocks linked to Scheme 2; spectra are source context."),
       ("heteroarylation general procedure", "Detailed common conditions/workup linked to the generic and eight Scheme 3 records."),
       ("4a-4c characterization", "Isolated masses for 4a-c linked to their scope experiments; analytical peaks are not reactions." )],
}
coverage = {"document_id":DOC, "semantic_accuracy_verified":False, "pages":[]}
unresolved_by_page = {
    1:["Graphical structures are generic Ar/R series and are not encoded as single concrete molecules."],
    2:["Prior-work Scheme 1 arrows omit yields and full variable identities; retained as generic literature examples."],
    3:["Gamma-substituted failure lacks substrate identity and conditions; no separate experiment fabricated.", "Absolute configuration of scope products is unreported."],
    4:["Absolute configurations are unreported; dr measurement methods/configurations for 6 and 7 are unreported.", "Mechanistic intermediates A-I are proposed, not isolated experiments."],
    5:["SI is not associated; cited precursor preparation is unavailable.", "Scheme 6 omits scale, solvent volume, workup/purification and yield basis."],
    6:["Characterization blocks do not establish full stereochemical composition.", "Heteroarylation products retain unassigned stereogenic centers."],
}
for p in range(1,7):
    evidence = [E[p], I[p]]
    coverage["pages"].append({
        "page":p, "state":"attempted",
        "regions":[{"region":region,"finding":finding,"evidence_ids":evidence} for region,finding in coverage_regions[p]],
        "record_ids":sorted(set(record_ids_by_page[p])),
        "unresolved":unresolved_by_page[p],
    })
(RUN / "coverage/pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

stereo_diags = {}
for mid in ["1a", "2a", "2f", "2r", "indole", "n_methylindole"]:
    smi = st.value(molecules[mid]["smiles"])
    if smi:
        stereo_diags[mid] = st.stereo_inspect(smi)
(RUN / "raw_model/stereo_diagnostics_pages_0001_0006.json").write_text(json.dumps(stereo_diags, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({"molecules":len(bundle["molecules"]), "reactions":len(bundle["reactions"]),
                  "coverage_pages":len(coverage["pages"])}, ensure_ascii=False))
