#!/usr/bin/env python3
"""Author the page 1-6 extraction bundle for Wang et al., ACS Catal. 2021.

The constants below are a transcription of the supplied page images.  This
script does not call a model, alter the schema, or certify chemical accuracy.
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
DOC = "doc_a3d4333b7270ebfb"
TITLE = "Enantioselective Construction of 1H-Isoindoles Containing Tri- and Difluoromethylated Quaternary Stereogenic Centers via Palladium-Catalyzed C-H Bond Imidoylation"
E = {p: [f"{DOC}:page{p}:image", f"{DOC}:page{p}:text"] for p in range(1, 7)}


def dump(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


def f(raw=None, normalized=None, evidence=(), method="original_extraction", status=None,
      verification=None, derivation=None, alternatives=()):
    if normalized is None and raw is not None and method == "original_extraction":
        normalized = raw
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   verification=verification, derivation=derivation,
                   alternatives=alternatives)


def q(number, unit, **extra):
    return {"number": number, "unit": unit, **extra}


bundle = json.loads((RUN / "submitted.json").read_text())
# The mechanical PDF preparation did not populate the title.  Add the title
# transcribed from page 1 to the prepared document metadata before merge, and
# keep any already-merged standard document copy synchronized.
dataset_path = RUN / "dataset.json"
dataset = json.loads(dataset_path.read_text())
for doc in dataset["documents"]:
    if doc["document_id"] == DOC:
        doc["title"] = TITLE
for doc in dataset.get("standard", {}).get("documents", []):
    if doc["document_id"] == DOC:
        doc["title"] = TITLE
dump(dataset_path, dataset)
for doc in bundle["documents"]:
    if doc["document_id"] == DOC:
        doc["title"] = TITLE
molecules = {}


def molecule(mid, label, evidence, *, name=None, smiles=None, stereo="unresolved",
             note=None, geometry=None):
    if mid in molecules:
        return mid
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["occurrences"][0]["evidence_ids"] = list(evidence)
    m["names"]["original"] = f(name or label, name or label, evidence)
    m["names"]["english"] = f(name or label, name or label, evidence)
    m["identity_note"] = f(note, note, evidence) if note else f()
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = f(
            f"candidate connectivity derived from the source drawing for {label}",
            smiles, evidence, method="structure_derivation",
            derivation="Connectivity was manually transcribed from the supplied source drawing; RDKit parsing is a graph diagnostic, not source verification.")
        if stereo in {"achiral", "not_applicable"}:
            m["stereochemistry"]["chirality_status"] = f(
                "no stereogenic element in the transcribed graph", "achiral", evidence,
                method="structure_derivation", derivation="RDKit found no potential atom or bond stereogenic element in the candidate graph; the source drawing was also inspected.")
            m["stereochemistry"]["geometry_status"] = f(
                "no source-relevant E/Z or ring cis/trans assignment", "not_applicable", evidence,
                method="structure_derivation", derivation="The candidate graph has no source-relevant configurational geometry.")
            m["stereochemistry"]["representation_status"] = f(
                "stereochemistry not applicable", "not_applicable", evidence,
                method="structure_derivation", derivation="No stereogenic feature needs encoding for this candidate graph.")
            m["stereochemistry"]["isomeric_smiles"] = f(
                None, diag["canonical_isomeric_smiles"], evidence, method="structure_derivation",
                derivation="RDKit canonical isomeric SMILES for the manually transcribed connectivity.")
            m["stereochemistry"]["identity_key"] = f(
                None, diag["stereo_graph_key"], evidence, method="structure_derivation",
                derivation="RDKit stereo graph key; not a sample identity or source verification.")
        elif stereo == "racemic":
            m["stereochemistry"]["atom_centers"] = f(
                "quaternary stereogenic carbon shown without an absolute descriptor",
                [{"center": "quaternary carbon bearing Ph and CF3", "configuration": "unspecified"}], evidence)
            m["stereochemistry"]["mixture"] = f("racemic product", {"type": "racemate", "ratio": "1:1"}, evidence)
            m["stereochemistry"]["chirality_status"] = f("racemic", "racemic", evidence)
            m["stereochemistry"]["geometry_status"] = f("not applicable", "not_applicable", evidence)
            m["stereochemistry"]["representation_status"] = f("stereochemistry not encoded in connectivity SMILES", "not_encoded", evidence)
        elif stereo in {"ee_mixture", "de_mixture"}:
            mixture_type = "enantiomeric_mixture" if stereo == "ee_mixture" else "diastereomeric_mixture"
            m["stereochemistry"]["atom_centers"] = f(
                "stereogenic center(s) are drawn but no R/S descriptor is printed",
                [{"configuration": "major configuration not assigned in pages 1-6"}], evidence)
            m["stereochemistry"]["mixture"] = f(
                "composition is experiment-specific and recorded as ee/de outcomes",
                {"type": mixture_type, "composition": "experiment-specific"}, evidence)
            m["stereochemistry"]["chirality_status"] = f(mixture_type, "mixture", evidence)
            m["stereochemistry"]["geometry_status"] = f(geometry or "not applicable", geometry or "not_applicable", evidence)
            m["stereochemistry"]["representation_status"] = f("source stereo is not encoded in the connectivity candidate", "not_encoded", evidence)
    else:
        m["representation_type"] = "unresolved"
        m["stereochemistry"]["chirality_status"] = f(
            "source drawing inspected; complete structure/stereo not encoded", stereo if stereo in {
                "specified", "partially_specified", "unspecified", "racemic", "mixture", "achiral", "unresolved", "not_applicable"
            } else "unresolved", evidence)
        m["stereochemistry"]["geometry_status"] = f(
            "source geometry not fully encoded", geometry or "unresolved", evidence)
        m["stereochemistry"]["representation_status"] = f(
            "complete structure not encoded", "not_encoded" if stereo in {"mixture", "racemic", "specified", "partially_specified", "unspecified"} else "unresolved", evidence)
        if stereo == "mixture":
            m["stereochemistry"]["atom_centers"] = f(
                "product drawing contains stereogenic center(s); no complete absolute assignment transcribed",
                [{"configuration": "unassigned"}], evidence)
            m["stereochemistry"]["mixture"] = f(
                "selectivity composition is experiment-specific", {"type": "stereoisomeric_mixture", "composition": "experiment-specific"}, evidence)
    molecules[mid] = m
    return mid


# Shared concrete graph candidates.  Product 3a uses connectivity only: the
# source-reported enantiomeric composition remains an outcome of each experiment.
M1A = molecule("mol_1a_cf3_isocyanide", "1a", E[2],
               name="(2,2,2-trifluoro-1-isocyanoethane-1,1-diyl)dibenzene",
               smiles="[C-]#[N+]C(c1ccccc1)(c1ccccc1)C(F)(F)F", stereo="achiral")
M1H = molecule("mol_1_cf2h_isocyanide", "1 (R = H)", E[3] + E[4],
               name="(2,2-difluoro-1-isocyanoethane-1,1-diyl)dibenzene",
               smiles="[C-]#[N+]C(c1ccccc1)(c1ccccc1)C(F)F", stereo="achiral")
M9 = molecule("mol_9_methyl_isocyanide", "9", E[5],
              name="methyl-substituted diaryl isocyanide 9",
              smiles="[C-]#[N+]C(C)(c1ccccc1)c1ccccc1", stereo="achiral")
PHI = molecule("mol_iodobenzene_2a", "2a", E[2] + E[3] + E[5], name="iodobenzene",
               smiles="Ic1ccccc1", stereo="achiral")
PHBR = molecule("mol_bromobenzene", "bromobenzene", E[2], name="bromobenzene",
                smiles="Brc1ccccc1", stereo="achiral")
P3A = molecule("mol_3a_enantiomeric_mixture", "3a", E[2] + E[3] + E[4],
               name="trifluoromethylated 1H-isoindole 3a",
               smiles="N1=C(c2ccccc2)c2ccccc2C1(c3ccccc3)C(F)(F)F", stereo="ee_mixture",
               note="Connectivity candidate only. The source reports experiment-specific ee values; the major absolute configuration is not assigned in the processed pages.")
P3A_RAC = molecule("mol_3a_racemate", "rac-3a", E[2], name="racemic 3a",
                   smiles="N1=C(c2ccccc2)c2ccccc2C1(c3ccccc3)C(F)(F)F", stereo="racemic")


def simple(mid, label, evidence, name=None):
    return molecule(mid, label, evidence, name=name or label, stereo="unresolved")


PDOAC = simple("mol_palladium_acetate", "Pd(OAc)2", E[2] + E[3] + E[4] + E[5], "palladium(II) acetate")
PDCL = simple("mol_palladium_chloride", "PdCl2", E[2], "palladium(II) chloride")
PDDBA = simple("mol_pd2dba3", "Pd2(dba)3", E[2], "tris(dibenzylideneacetone)dipalladium(0)")
PPH3 = simple("mol_triphenylphosphine", "PPh3", E[2], "triphenylphosphine")
L1 = simple("mol_ligand_l1", "L1", E[2], "phosphoramidite ligand L1")
L2 = simple("mol_ligand_l2", "L2", E[2], "phosphoramidite ligand L2")
L3 = simple("mol_ligand_l3", "L3", E[2] + E[3] + E[4] + E[5], "SPINOL-derived phosphoramidite ligand L3")
BINAP = simple("mol_binap", "BINAP", E[2], "BINAP")
CS2CO3 = simple("mol_cs2co3", "Cs2CO3", E[2] + E[3] + E[4] + E[5], "cesium carbonate")
PIVOH = simple("mol_pivoh", "PivOH", E[2] + E[3] + E[4] + E[5], "pivalic acid")
TOL = simple("mol_toluene", "toluene", E[2] + E[3] + E[4] + E[5], "toluene")
HNO3 = simple("mol_hno3", "HNO3", E[4], "nitric acid")
H2SO4 = simple("mol_h2so4", "H2SO4", E[4], "sulfuric acid")


def product(mid, label, evidence, mode="ee", note=None):
    stereo = "mixture"
    return molecule(mid, label, evidence, name=f"source-identified product {label}", stereo=stereo,
                    geometry="unresolved" if label in {"3y", "3z", "4x", "4y"} else "not_applicable",
                    note=note or "Complete connectivity and absolute configuration were not encoded; retain the source drawing and experiment-specific selectivity outcome.")


def partner_for(label, evidence, kind="iodide"):
    return simple(f"mol_partner_{label}", f"{kind} coupling partner for {label}", evidence,
                  f"source-depicted {kind} coupling partner leading to {label}")


def modified_isocyanide_for(label, evidence):
    return simple(f"mol_isocyanide_{label}", f"isocyanide precursor for {label}", evidence,
                  f"source-depicted substituted isocyanide precursor leading to {label}")


def set_charge(p, key, raw, norm, evidence, method="original_extraction", derivation=None):
    p["charges"][key] = f(raw, norm, evidence, method=method, derivation=derivation)


def participant(pid, mid, role, evidence, *, amount=None, equivalents=None, mol_percent=None,
                volume=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, evidence)
    if amount is not None:
        set_charge(p, "amount", f"{amount} mmol", q(amount, "mmol"), evidence)
    if equivalents is not None:
        set_charge(p, "equivalents", f"{equivalents} equiv", q(equivalents, "equiv"), evidence)
    if mol_percent is not None:
        set_charge(p, "mol_percent", f"{mol_percent} mol %", q(mol_percent, "mol%"), evidence)
    if volume is not None:
        set_charge(p, "volume", volume[0], q(volume[1], "mL"), evidence)
    if reference is not None:
        p["equivalence_reference"] = f(reference, reference, evidence)
    if product_role:
        p["product_role"] = f(product_role, product_role, evidence)
    return p


def base_reaction(rid, evidence, category, location, name="palladium-catalyzed enantioselective C-H bond imidoylation"):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = f(category, [category], evidence)
    r["location"] = f(location, {"pages": sorted({int(e.split(':page')[1].split(':')[0]) for e in evidence if ':page' in e})}, evidence)
    r["reaction_name"] = f(name, name, evidence)
    r["transformation_type"] = f(
        None, "imidoylative cyclization to a 1H-isoindole", evidence,
        method="inference", derivation="Transformation type was inferred from the source reaction title, arrows, and product scaffold.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["No supporting-information file was supplied or explicitly associated; detailed workup, purification, and characterization remain unavailable.",
                                  "Complete atom mapping and absolute stereochemical descriptors remain unresolved."]}
    r["completeness"] = "partial"
    return r


def add_step(r, evidence, participant_ids, *, operation, temperature, time, atmosphere="Ar"):
    s = st.new_step("reaction", evidence)
    s["participant_ids"] = participant_ids
    s["operation"] = f(operation, operation, evidence)
    s["temperature"] = f(temperature if isinstance(temperature, str) else f"{temperature} °C",
                           temperature if isinstance(temperature, str) else q(temperature, "°C"), evidence)
    s["time"] = f(time if isinstance(time, str) else f"{time} h",
                    time if isinstance(time, str) else q(time, "h"), evidence)
    s["atmosphere"] = f(atmosphere, atmosphere, evidence) if atmosphere else f()
    r["steps"] = [s]


def classify_imidoylation(r, evidence, c_c_count=2):
    changes = []
    for i in range(c_c_count):
        changes.append({"change_id": f"b{i+1}", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
                        "before_order": 0, "after_order": 1, "change_type": "formed"})
    j = len(changes) + 1
    changes.extend([
        {"change_id": f"b{j}", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 3, "after_order": 2, "change_type": "order_decreased"},
        {"change_id": f"b{j+1}", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "I"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
        {"change_id": f"b{j+2}", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "hydrogen_count": 1},
    ])
    r["classification"]["bond_changes"] = [f(None, ch, evidence, method="structure_derivation",
                                                      derivation="Partial net-bond interpretation from the source substrate/product drawings; atom correspondence is not fully mapped.") for ch in changes]
    r["classification"]["labels"] = f(None, [{"category_number": 2, "code": "C-C", "supporting_change_ids": [f"b{i+1}" for i in range(c_c_count)]}], evidence,
                                             method="structure_derivation", derivation="Only newly formed C-C connections receive a formation label; C-N order decrease and C-I/C-H losses are retained but not labeled as new bonds.")
    r["classification"]["coverage"] = f(None, "partial", evidence, method="inference",
                                               derivation="Complete atom mapping, proton balance, and coproduct identities are unavailable in the processed pages.")
    r["classification"]["limitations"] = ["Partial net transformation only; no complete atom mapping.", "Implicit hydrogen loss and halide departure are not assigned to specific mapped atoms."]


def outcome(r, oid, kind, target, raw, norm, evidence, *, isolated=False, metric=None, scope_extra=None):
    o = st.new_outcome(oid, kind, [target], evidence)
    o["value"] = f(raw, norm, evidence)
    if kind == "yield":
        if isolated:
            o["yield_type"] = f("isolated yield", "isolated", evidence)
            o["measurement_basis"] = f("Yields of isolated products are given.", "isolated", evidence)
        else:
            o["yield_type"] = f()
            o["measurement_basis"] = f()
    elif kind == "selectivity":
        scope = {"metric": metric, "major_configuration": "unassigned"}
        if scope_extra:
            scope.update(scope_extra)
        o["scope"] = f(metric, scope, evidence)
        o["measurement_basis"] = f()
    elif kind == "product_amount":
        o["measurement_basis"] = f("reported product mass", "reported mass", evidence)
    r["outcomes"].append(o)


def standard_components(r, evidence, substrate, partner, product_mid, *, scale, partner_amount,
                        pd=PDOAC, pd_load=1, ligand=L3, ligand_load=2, cs_amount=None,
                        cs_equiv=None, piv_amount=None, piv_equiv=None, toluene_ml=None,
                        product_role="main_product"):
    ref = substrate
    ps = [participant("substrate", substrate, "reactant", evidence, amount=scale),
          participant("partner", partner, "reactant", evidence, amount=partner_amount),
          participant("pd", pd, "catalyst", evidence, mol_percent=pd_load, reference=ref),
          participant("ligand", ligand, "ligand", evidence, mol_percent=ligand_load, reference=ref)]
    if cs_amount is not None:
        ps.append(participant("cs2co3", CS2CO3, "reagent", evidence, amount=cs_amount))
    elif cs_equiv is not None:
        ps.append(participant("cs2co3", CS2CO3, "reagent", evidence, equivalents=cs_equiv, reference=ref))
    if piv_amount is not None:
        ps.append(participant("pivoh", PIVOH, "additive", evidence, amount=piv_amount))
    elif piv_equiv is not None:
        ps.append(participant("pivoh", PIVOH, "additive", evidence, equivalents=piv_equiv, reference=ref))
    if toluene_ml is not None:
        ps.append(participant("toluene", TOL, "solvent", evidence, volume=(f"{toluene_ml} mL total", toluene_ml)))
        r["solvent_conditions"]["composition"] = f("toluene", ["toluene"], evidence)
        r["solvent_conditions"]["total_volume"] = f(f"{toluene_ml} mL", q(toluene_ml, "mL"), evidence)
    else:
        ps.append(participant("toluene", TOL, "solvent", evidence))
        r["solvent_conditions"]["composition"] = f("toluene", ["toluene"], evidence)
    ps.append(participant("product", product_mid, "product", evidence, product_role=product_role))
    r["participants"] = ps
    r["reaction_scale"] = f(f"{scale} mmol {substrate}", {"reference_molecule_id": substrate, "amount": q(scale, "mmol")}, evidence)
    add_step(r, evidence, [p["participant_id"] for p in ps if p["stage"] == "reaction"],
             operation="isocyanide solution added dropwise by syringe pump over the reported reaction period",
             temperature=60, time=2, atmosphere="Ar")
    classify_imidoylation(r, evidence)


reactions = []

# Figure 1b: cited prior work, and Figure 1c/page-1 graphical abstract: this-work generic scheme.
lit_sub = simple("mol_generic_cramer_imidoyl_chloride", "Cramer imidoyl chloride substrate", E[2])
lit_prod = simple("mol_generic_cramer_isoindole", "Cramer isoindole product", E[2])
lit = st.new_reaction(f"{DOC}:figure1b_cramer_literature", DOC, E[2])
lit["record_type"] = "literature_example"
lit["experiment_categories"] = f("cited previous work", ["other"], E[2])
lit["location"] = f("Figure 1b", {"page": 2, "figure": "1b"}, E[2])
lit["reaction_name"] = f("Cramer's intramolecular synthesis of 1H-isoindoles", "Cramer's intramolecular synthesis of 1H-isoindoles", E[2])
lit["transformation_type"] = f("literature example", "intramolecular C-H bond imidoylation", E[2])
lit["participants"] = [participant("substrate", lit_sub, "reactant", E[2]), participant("product", lit_prod, "product", E[2], product_role="main_product")]
lit["classification"]["coverage"] = f(None, "unresolved", E[2], method="inference", derivation="The cited generic scheme was not atom mapped in this extraction.")
lit["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [], "unresolved": ["Cited prior work; generic R/Ar/RF structures were not expanded."]}
lit["completeness"] = "partial"
reactions.append(lit)

gen_sub1 = simple("mol_generic_isocyanide_1", "generic isocyanide 1", E[1] + E[2])
gen_sub2 = simple("mol_generic_iodide_2", "generic aryl/alkenyl iodide 2", E[1] + E[2])
gen_prod = product("mol_generic_isoindole_3_4", "generic product 3/4", E[1] + E[2])
gen = st.new_reaction(f"{DOC}:figure1c_current_generic", DOC, E[1] + E[2])
gen["record_type"] = "generic_scheme"
gen["experiment_categories"] = f("this work generic scheme", ["other"], E[1] + E[2])
gen["location"] = f("page 1 graphical abstract and Figure 1c", {"pages": [1, 2]}, E[1] + E[2])
gen["reaction_name"] = f("Intermolecular access to 1H-isoindoles with CF3/CF2H-containing quaternary stereogenic centers", "Intermolecular access to 1H-isoindoles", E[1] + E[2])
gen["transformation_type"] = f(None, "generic imidoylative cyclization", E[1] + E[2], method="inference", derivation="Inferred from the generic source arrow.")
gen["participants"] = [participant("isocyanide", gen_sub1, "reactant", E[1] + E[2]), participant("iodide", gen_sub2, "reactant", E[1] + E[2]), participant("product", gen_prod, "product", E[1] + E[2], product_role="main_product")]
classify_imidoylation(gen, E[1] + E[2])
gen["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [], "unresolved": ["Generic Ar/RF/tether variables and ligand L* are not expanded."]}
gen["completeness"] = "partial"
reactions.append(gen)
GENERAL = gen["reaction_id"]


# Table 1 optimization. Entry 9 is the same 0.5 mmol standard-scope 3a experiment
# repeated in Scheme 1, so it is authored once below with both evidence sources.
table_rows = [
    (1, PDOAC, PPH3, 10, 20, 97, None, P3A_RAC),
    (2, PDOAC, L1, 10, 20, 96, 92, P3A),
    (3, PDOAC, L2, 10, 20, 91, 11, P3A),
    (4, PDOAC, L3, 10, 20, 94, 98, P3A),
    (5, PDOAC, BINAP, 10, 20, 30, 8, P3A),
    (6, PDCL, L3, 10, 20, 60, 95, P3A),
    (7, PDDBA, L3, 10, 20, 49, 99, P3A),
    (8, PDOAC, L3, 5, 10, 93, 97, P3A),
]
for entry, pd, lig, pdl, ligl, yld, ee, target in table_rows:
    ev = E[2]
    r = base_reaction(f"table1_entry_{entry:02d}", ev, "optimization", f"Table 1, entry {entry}")
    r["general_scheme_ids"] = [GENERAL]
    standard_components(r, ev, M1A, PHI, target, scale=0.1, partner_amount=0.12, pd=pd, pd_load=pdl,
                        ligand=lig, ligand_load=ligl, cs_amount=0.1, piv_amount=0.06, toluene_ml=2.0)
    outcome(r, "yield", "yield", target, f"{yld}%", q(yld, "%"), ev, isolated=True)
    if ee is not None:
        outcome(r, "ee", "selectivity", target, f"{ee}% ee", q(ee, "%"), ev, metric="ee")
    else:
        r["notes"].append(f("racemic product 3a", "racemic product 3a", ev))
    reactions.append(r)

# Table 1 entry 9 / Scheme 1 product 3a.  The duplicated presentation has
# conflicting/ambiguous absolute charges and solvent volumes, retained below.
ev_3a = E[2] + E[3]
r = base_reaction("table1_entry_09_scope_3a", ev_3a, "optimization", "Table 1 entry 9 and Scheme 1 product 3a")
r["experiment_categories"] = f("optimization and substrate scope", ["optimization", "substrate_scope"], ev_3a)
r["general_scheme_ids"] = [GENERAL]
standard_components(r, ev_3a, M1A, PHI, P3A, scale=0.5, partner_amount=0.75, pd_load=1, ligand_load=2,
                    cs_amount=0.5, piv_amount=0.3, toluene_ml=6.0)
def conflict_fact(a_raw, a_norm, a_ev, b_raw, b_norm, b_ev):
    return st.fact([a_raw, b_raw], None, list(dict.fromkeys(a_ev + b_ev)),
                   status="conflict",
                   alternatives=[f(a_raw, a_norm, a_ev), f(b_raw, b_norm, b_ev)])
r["participants"][1]["charges"]["amount"] = conflict_fact("0.12 mmol (Table 1 general footnote)", q(0.12, "mmol"), E[2], "0.75 mmol (Scheme 1)", q(0.75, "mmol"), E[3])
r["participants"][4]["charges"]["amount"] = conflict_fact("0.1 mmol (Table 1 general footnote)", q(0.1, "mmol"), E[2], "0.5 mmol (Scheme 1)", q(0.5, "mmol"), E[3])
r["participants"][5]["charges"]["amount"] = conflict_fact("0.06 mmol (Table 1 general footnote)", q(0.06, "mmol"), E[2], "0.3 mmol (Scheme 1)", q(0.3, "mmol"), E[3])
r["participants"][6]["charges"]["volume"] = conflict_fact("1.0 mL + 1.0 mL (Table 1 general footnote)", q(2.0, "mL"), E[2], "3 mL + 3 mL (Scheme 1)", q(6.0, "mL"), E[3])
r["solvent_conditions"]["total_volume"] = conflict_fact("2.0 mL (Table 1 general footnote)", q(2.0, "mL"), E[2], "6 mL (Scheme 1)", q(6.0, "mL"), E[3])
r["notes"].append(f("The repeated 93% yield/97% ee experiment has inconsistent or underspecified absolute charges between Table 1 footnote c and Scheme 1; both presentations are retained.",
                    "source presentation conflict retained", ev_3a))
r["review"]["unresolved"].append("Resolve whether Table 1 footnote c rescales all general charges/volumes; the printed text only explicitly changes reaction scale and catalyst/ligand loading.")
outcome(r, "yield", "yield", P3A, "93%", q(93, "%"), ev_3a, isolated=True)
outcome(r, "ee", "selectivity", P3A, "97% ee", q(97, "%"), ev_3a, metric="ee")
reactions.append(r)

# Table 1 entry 10: bromobenzene electrophile, 0.5 mmol scale footnote c.
r = base_reaction("table1_entry_10", E[2], "optimization", "Table 1, entry 10")
r["general_scheme_ids"] = [GENERAL]
standard_components(r, E[2], M1A, PHBR, P3A, scale=0.5, partner_amount=0.12, pd_load=1, ligand_load=2,
                    cs_amount=0.1, piv_amount=0.06, toluene_ml=2.0)
r["notes"].append(f("Footnote c states a 0.5 mmol scale but does not explicitly restate scaled bromobenzene/additive/solvent charges; printed general-footnote values are retained without recalculation.",
                    "scale/charge relationship unresolved", E[2]))
r["review"]["unresolved"].append("Entry 10 scale is explicit, but scaling of the general-footnote absolute charges is not explicit.")
outcome(r, "yield", "yield", P3A, "92%", q(92, "%"), E[2], isolated=True)
outcome(r, "ee", "selectivity", P3A, "89% ee", q(89, "%"), E[2], metric="ee")
reactions.append(r)


# Scheme 1 CF3 scope. The standard 3a row above is not duplicated here.
cf3 = [
    ("3b",90,89,None),("3c",94,96,None),("3d",91,97,None),("3e",98,96,None),
    ("3f",93,90,None),("3g",85,93,None),("3h",95,95,None),
    ("3i",96,99,"b"),("3j",94,98,None),("3k",88,93,None),
    ("3l",95,97,None),("3m",91,97,None),("3n",97,98,None),("3o",94,96,None),
    ("3p",99,95,None),("3q",97,92,None),("3r",90,97,None),
    ("3s",91,97,"bc"),("3t",93,99,"bc"),("3u",78,99,None),
    ("3v",93,92,None),("3w",95,94,None),("3x",81,98,None),
    ("3y",82,99,"b"),("3z",57,99,"b"),
    ("3aa",99,98,"d"),("3ab",80,96,"d"),("3ac",94,98,"d"),("3ad",92,96,"d"),("3ae",82,96,"d"),
]
for label, yld, ee, flag in cf3:
    ev = E[3]
    prod = product(f"mol_{label}", label, ev)
    if label in {"3aa", "3ab", "3ac", "3ad", "3ae"}:
        sub, partner = modified_isocyanide_for(label, ev), PHI
    else:
        sub, partner = M1A, partner_for(label, ev)
    pd_load, lig_load, temp = 1, 2, 60
    if flag and "b" in flag: pd_load, lig_load = 10, 20
    if flag and "d" in flag: pd_load, lig_load = 5, 10
    if flag and "c" in flag: temp = 70
    r = base_reaction(f"scheme1_scope_{label}", ev, "substrate_scope", f"Scheme 1, product {label}")
    r["general_scheme_ids"] = [GENERAL]
    standard_components(r, ev, sub, partner, prod, scale=0.5, partner_amount=0.75, pd_load=pd_load,
                        ligand_load=lig_load, cs_amount=0.5, piv_amount=0.3, toluene_ml=6.0)
    r["steps"][0]["temperature"] = f(f"{temp} °C", q(temp, "°C"), ev)
    outcome(r, "yield", "yield", prod, f"{yld}%", q(yld, "%"), ev, isolated=True)
    outcome(r, "ee", "selectivity", prod, f"{ee}% ee", q(ee, "%"), ev, metric="ee")
    reactions.append(r)

# Two additional 3a scale demonstrations displayed separately below product 3a.
for suffix, raw_scale, yld, ee, mass in [
    ("5mmol", 5.0, 92, None, None),
    ("gram", None, 99, 90, 1.67),
]:
    r = base_reaction(f"scheme1_scope_3a_{suffix}", E[3], "scale_up", f"Scheme 1, 3a {suffix} scale demonstration")
    r["general_scheme_ids"] = [GENERAL]
    # Only explicitly printed scale/product mass is structured; no proportional charges are calculated.
    ps = [participant("substrate", M1A, "reactant", E[3], amount=raw_scale) if raw_scale else participant("substrate", M1A, "reactant", E[3]),
          participant("partner", PHI, "reactant", E[3]),
          participant("pd", PDOAC, "catalyst", E[3], mol_percent=1, reference=M1A),
          participant("ligand", L3, "ligand", E[3], mol_percent=2, reference=M1A),
          participant("cs2co3", CS2CO3, "reagent", E[3], equivalents=1, reference=M1A),
          participant("pivoh", PIVOH, "additive", E[3], equivalents=0.6, reference=M1A),
          participant("toluene", TOL, "solvent", E[3]),
          participant("product", P3A, "product", E[3], product_role="main_product")]
    r["participants"] = ps
    if raw_scale:
        r["reaction_scale"] = f("5 mmol scale", {"reference_molecule_id": M1A, "amount": q(5, "mmol")}, E[3])
    r["solvent_conditions"]["composition"] = f("toluene", ["toluene"], E[3])
    add_step(r, E[3], [p["participant_id"] for p in ps], operation="conditions inherited from Scheme 1; isocyanide added by syringe pump", temperature=60, time=2)
    classify_imidoylation(r, E[3])
    outcome(r, "yield", "yield", P3A, f"{yld}%", q(yld, "%"), E[3], isolated=True)
    if ee is not None:
        outcome(r, "ee", "selectivity", P3A, f"{ee}% ee", q(ee, "%"), E[3], metric="ee")
    if mass is not None:
        outcome(r, "product_mass", "product_amount", P3A, f"{mass} g", q(mass, "g"), E[3])
    r["review"]["unresolved"].append("The scale-up callout does not restate all absolute reagent and solvent charges; no proportional amounts were calculated.")
    reactions.append(r)


# Scheme 1 CF2H scope, including the explicitly discussed low-catalyst 4a result.
P4A = product("mol_4a", "4a", E[3] + E[4], note="The source includes an X-ray model (CCDC 2086996), but no R/S descriptor was transcribed from pages 1-6; the ee is experiment-specific.")
r = base_reaction("scheme1_scope_4a_low_catalyst", E[3] + E[4], "optimization", "Scheme 1/prose, 4a under standard 1 mol% catalyst conditions")
r["general_scheme_ids"] = [GENERAL]
standard_components(r, E[3] + E[4], M1H, PHI, P4A, scale=0.5, partner_amount=0.75, pd_load=1, ligand_load=2,
                    cs_amount=0.5, piv_amount=0.3, toluene_ml=6.0)
outcome(r, "yield", "yield", P4A, "43%", q(43, "%"), E[3] + E[4], isolated=True)
outcome(r, "ee", "selectivity", P4A, "87% ee", q(87, "%"), E[3] + E[4], metric="ee")
r["notes"].append(f("a large amount of isocyanide remained unreacted", "large amount of isocyanide remaining", E[4]))
reactions.append(r)

cf2 = [
    ("4a",95,92),("4b",86,95),("4c",85,96),("4d",93,92),("4e",89,89),
    ("4f",92,94),("4g",96,92),("4h",89,96),("4i",95,97),("4j",93,96),("4k",88,94),
    ("4l",90,97),("4m",95,91),("4n",91,86),("4o",90,96),("4p",95,97),
    ("4q",97,96),("4r",97,97),("4s",97,97),("4t",93,97),("4u",91,95),
    ("4v",87,95),("4w",95,94),("4x",77,97),("4y",24,96),
    ("4z",91,93),("4aa",84,98),("4ab",80,96),
]
for label, yld, ee in cf2:
    ev = E[3] + (E[4] if label in {"4a", "4x", "4y", "4z", "4aa", "4ab"} else [])
    prod = P4A if label == "4a" else product(f"mol_{label}", label, ev)
    if label in {"4z", "4aa", "4ab"}:
        sub, partner = modified_isocyanide_for(label, ev), PHI
    else:
        sub, partner = M1H, partner_for(label, ev)
    r = base_reaction(f"scheme1_scope_{label}", ev, "substrate_scope", f"Scheme 1, product {label}")
    r["general_scheme_ids"] = [GENERAL]
    standard_components(r, ev, sub, partner, prod, scale=0.1, partner_amount=0.15, pd_load=5, ligand_load=10,
                        cs_amount=0.1, piv_amount=0.06, toluene_ml=None)
    outcome(r, "yield", "yield", prod, f"{yld}%", q(yld, "%"), ev, isolated=True)
    outcome(r, "ee", "selectivity", prod, f"{ee}% ee", q(ee, "%"), ev, metric="ee")
    if label == "4y":
        byp = product("mol_4y_prime", "4y'", ev, note="Major imidic pivalic anhydride byproduct; complete connectivity and stereochemistry were not encoded.")
        r["participants"].append(participant("byproduct", byp, "product", ev, product_role="major_byproduct"))
        outcome(r, "byproduct_yield", "yield", byp, "56%", q(56, "%"), ev, isolated=True)
        r["notes"].append(f("4y' was the major byproduct", "imidic pivalic anhydride 4y' major byproduct", ev))
        r["review"]["unresolved"].append("Bond changes for byproduct 4y' were not classified separately from the main-product transformation.")
    reactions.append(r)


# Scheme 2a allene insertion cascade: three distinct arrows/products.
for label, substrate_label, yld, ee, scaffold in [
    ("6a", "5a", 84, 98, "benzofuran-isoindole"),
    ("6b", "5b", 54, 98, "indole-isoindole"),
    ("6c", "5c", 77, 99, "isoquinolinone-isoindole"),
]:
    ev = E[4]
    sub5 = simple(f"mol_{substrate_label}", substrate_label, ev, f"source-depicted allene substrate {substrate_label}")
    prod = product(f"mol_{label}", label, ev, note=f"Enantioenriched {scaffold} product; source reports ee but no absolute descriptor in pages 1-6.")
    r = base_reaction(f"scheme2a_{label}", ev, "substrate_scope", f"Scheme 2a, {substrate_label} to {label}",
                      name="allene insertion/isocyanide insertion cascade")
    r["general_scheme_ids"] = [GENERAL]
    standard_components(r, ev, M1A, sub5, prod, scale=0.1, partner_amount=0.15, pd_load=10, ligand_load=20,
                        cs_amount=0.1, piv_amount=0.06, toluene_ml=None)
    r["transformation_type"] = f(None, f"allene insertion cascade to {scaffold}", ev, method="inference",
                                      derivation="Inferred from Scheme 2a and the authors' description of intramolecular allene insertion followed by isocyanide insertion.")
    outcome(r, "yield", "yield", prod, f"{yld}%", q(yld, "%"), ev, isolated=True)
    outcome(r, "ee", "selectivity", prod, f"{ee}% ee", q(ee, "%"), ev, metric="ee")
    r["classification"]["limitations"].append("The additional allene-cascade bond changes are not completely atom mapped; the C-C label is partial.")
    reactions.append(r)


# Scheme 2b natural-product/drug/dye-derived iodides. Pairs share the pictured iodide.
scheme2b = [
    ("7a", M1A, "7ab", 95, 93, "de"),("7b", M1H, "7ab", 99, 96, "de"),
    ("7c", M1A, "7cd", 96, 99, "de"),("7d", M1H, "7cd", 94, 97, "de"),
    ("7e", M1A, "7ef", 96, 97, "de"),("7f", M1H, "7ef", 85, 93, "de"),
    ("7g", M1A, "7g", 65, 85, "ee"),("7h", M1A, "7h", 72, 89, "de"),
]
partner_cache = {}
for label, sub, pkey, yld, sel, metric in scheme2b:
    ev = E[4]
    partner = partner_cache.setdefault(pkey, partner_for(pkey, ev))
    prod = product(f"mol_{label}", label, ev, note=f"Complex source-derived product {label}; {metric} is recorded per experiment, while full inherited/new stereochemistry and absolute assignment remain unencoded.")
    r = base_reaction(f"scheme2b_{label}", ev, "substrate_scope", f"Scheme 2b, product {label}")
    r["general_scheme_ids"] = [GENERAL]
    # Scheme 2b reports relative loadings but no explicit substrate mmol.
    ps = [participant("substrate", sub, "reactant", ev), participant("partner", partner, "reactant", ev),
          participant("pd", PDOAC, "catalyst", ev, mol_percent=5, reference=sub),
          participant("ligand", L3, "ligand", ev, mol_percent=10, reference=sub),
          participant("cs2co3", CS2CO3, "reagent", ev, equivalents=1, reference=sub),
          participant("pivoh", PIVOH, "additive", ev, equivalents=0.6, reference=sub),
          participant("toluene", TOL, "solvent", ev), participant("product", prod, "product", ev, product_role="main_product")]
    r["participants"] = ps
    r["solvent_conditions"]["composition"] = f("toluene", ["toluene"], ev)
    add_step(r, ev, [p["participant_id"] for p in ps], operation="reaction under the Scheme 2b conditions", temperature=60, time=2)
    classify_imidoylation(r, ev)
    outcome(r, "yield", "yield", prod, f"{yld}%", q(yld, "%"), ev, isolated=True)
    outcome(r, metric, "selectivity", prod, f"{sel}% {metric}", q(sel, "%"), ev, metric=metric)
    reactions.append(r)


# Scheme 2c mono- and dinitration of 3a are separate experiments.
for label, yld, acid_v, sulfur_v, hours, nitro_count in [
    ("8a", 56, 0.4, 0.26, 16, 1),
    ("8b", 52, 2.0, 1.3, 32, 2),
]:
    ev = E[4]
    prod = product(f"mol_{label}", label, ev, note=f"Nitrated derivative {label}; the source reports retention of 97% ee, but the absolute configuration is not encoded.")
    r = base_reaction(f"scheme2c_{label}", ev, "derivatization", f"Scheme 2c, 3a to {label}", name="selective nitration of 3a")
    r["transformation_type"] = f("mono- or dinitration", "aromatic nitration", ev)
    ps = [participant("substrate", P3A, "reactant", ev, amount=0.2),
          participant("hno3", HNO3, "reagent", ev, volume=(f"{acid_v} mL", acid_v)),
          participant("h2so4", H2SO4, "reagent", ev, volume=(f"{sulfur_v} mL", sulfur_v)),
          participant("product", prod, "product", ev, product_role="main_product")]
    r["participants"] = ps
    r["reaction_scale"] = f("3a (0.2 mmol)", {"reference_molecule_id": P3A, "amount": q(0.2, "mmol")}, ev)
    add_step(r, ev, [p["participant_id"] for p in ps], operation="nitration in mixed nitric and sulfuric acids", temperature="rt", time=hours, atmosphere=None)
    # Replace imidoylation classification with nitration C-N formations.
    changes = []
    for i in range(nitro_count):
        changes.append({"change_id": f"b{i+1}", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"})
        changes.append({"change_id": f"h{i+1}", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "hydrogen_count": 1})
    r["classification"]["bond_changes"] = [f(None, ch, ev, method="structure_derivation", derivation="Partial net aromatic nitration change inferred from Scheme 2c drawings.") for ch in changes]
    r["classification"]["labels"] = f(None, [{"category_number": 5, "code": "C-N", "supporting_change_ids": [f"b{i+1}" for i in range(nitro_count)]}], ev,
                                             method="structure_derivation", derivation="Each newly installed nitro group creates one aryl C-N connection; displaced C-H bonds are retained as losses.")
    r["classification"]["coverage"] = f(None, "partial", ev, method="inference", derivation="No complete atom mapping or acid-derived coproduct mapping is available.")
    r["classification"]["limitations"] = ["Partial net nitration only; nitro oxygen provenance and complete atom mapping unresolved."]
    outcome(r, "yield", "yield", prod, f"{yld}%", q(yld, "%"), ev, isolated=True)
    outcome(r, "ee", "selectivity", prod, "97% ee", q(97, "%"), ev, metric="ee")
    reactions.append(r)


# Scheme 3a true control experiment.
P10 = product("mol_10", "10", E[5], note="Control imidoylation product with a methyl group instead of fluoroalkyl; source reports 57% ee with no absolute assignment.")
r = base_reaction("scheme3a_control_10", E[5], "control", "Scheme 3a, control experiment")
r["general_scheme_ids"] = [GENERAL]
standard_components(r, E[5], M9, PHI, P10, scale=0.5, partner_amount=0.75, pd_load=1, ligand_load=2,
                    cs_amount=0.5, piv_amount=0.3, toluene_ml=None)
outcome(r, "yield", "yield", P10, "11%", q(11, "%"), E[5], isolated=False)
outcome(r, "ee", "selectivity", P10, "57% ee", q(57, "%"), E[5], metric="ee")
r["notes"].append(f("fluoroalkyl group might play a crucial role in enantioselectivity control and reaction efficiency", "author interpretation", E[5]))
r["review"]["unresolved"].append("The control-scheme yield measurement basis is not stated in the processed pages.")
reactions.append(r)


# Scheme 3b is an author-proposed stereocontrol model, not a concrete experiment.
mech = st.new_reaction(f"{DOC}:scheme3b_proposed_mode", DOC, E[5])
mech["record_type"] = "proposed_mechanism"
mech["experiment_categories"] = f("catalytic mode proposal", ["mechanistic_study"], E[5])
mech["location"] = f("Scheme 3b", {"page": 5, "scheme": "3b"}, E[5])
mech["reaction_name"] = f("Enantiocontrol mode proposal", "Enantiocontrol mode proposal", E[5])
mech["transformation_type"] = f("plausible C-H···F-C hydrogen-bond effect", "proposed stereocontrol model", E[5])
mech["classification"]["coverage"] = f(None, "not_applicable", E[5], method="inference", derivation="Scheme 3b compares proposed transition-state modes and is not a net reaction experiment.")
mech["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [], "unresolved": ["Author-proposed model; no experimental outcome is assigned to the depicted favored/disfavored modes."]}
mech["completeness"] = "partial"
reactions.append(mech)


# Source context separates non-synthetic observations and document-level facts.
bundle["source_context"] = [
    {"context_id": f"{DOC}:article_metadata", "document_id": DOC, "kind": "bibliographic_metadata",
     "title": TITLE,
     "doi": "10.1021/acscatal.1c03682", "evidence_ids": E[1]},
    {"context_id": f"{DOC}:figure1a_bioactive_examples", "document_id": DOC, "kind": "background_figure",
     "description": "Figure 1a shows selected fluorinated drugs and bioactive molecules; it is not an experiment performed in this paper.", "evidence_ids": E[2]},
    {"context_id": f"{DOC}:xray_models", "document_id": DOC, "kind": "characterization",
     "description": "X-ray renderings are shown for 3u (CCDC 2019847), 4a (CCDC 2086996), 8a (CCDC 2087012), and 8b (CCDC 2086999). They support source characterization but are not separate synthesis experiments.", "evidence_ids": E[3] + E[4]},
    {"context_id": f"{DOC}:si_statement", "document_id": DOC, "kind": "supplement_statement",
     "description": "Page 5 states that Supporting Information contains experimental, characterization, and X-ray data, but no SI file was supplied or explicitly associated in this workspace.", "evidence_ids": E[5]},
    {"context_id": f"{DOC}:page6_references", "document_id": DOC, "kind": "bibliography",
     "description": "Page 6 contains references only; cited papers were not converted into current-paper experiments.", "evidence_ids": E[6]},
]

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
dump(RUN / "submitted.json", bundle)


# Raw transcription and page coverage are deliberately separate from the formal bundle.
raw = {
    "document_id": DOC,
    "pages_inspected": [1, 2, 3, 4, 5, 6],
    "chemical_accuracy_verified": False,
    "table1_rows": [{"entry": row[0], "yield_percent": row[5], "ee_percent": row[6]} for row in table_rows] +
                   [{"entry": 9, "yield_percent": 93, "ee_percent": 97}, {"entry": 10, "yield_percent": 92, "ee_percent": 89}],
    "scheme1_cf3": [{"label": a, "yield_percent": b, "ee_percent": c, "footnote_flags": d} for a,b,c,d in cf3],
    "scheme1_cf2h": [{"label": a, "yield_percent": b, "ee_percent": c} for a,b,c in cf2],
    "scheme2a": [{"product": a, "substrate": b, "yield_percent": c, "ee_percent": d} for a,b,c,d,_ in [("6a","5a",84,98,""),("6b","5b",54,98,""),("6c","5c",77,99,"")]],
    "scheme2b": [{"label": a, "yield_percent": d, "selectivity_percent": e, "metric": g} for a,_,_,d,e,g in scheme2b],
    "scheme2c": [{"label": "8a", "yield_percent": 56, "ee_percent": 97}, {"label": "8b", "yield_percent": 52, "ee_percent": 97}],
    "scheme3a": {"label": "10", "yield_percent": 11, "ee_percent": 57},
    "deduplication": ["Table 1 entry 9 and Scheme 1 3a (93% yield, 97% ee) are represented once with both page-2 and page-3 evidence."],
    "limitations": ["No SI was supplied.", "Many scope structures and absolute configurations remain source-drawing identities rather than complete SMILES.", "Machine-authored records remain needs_review."],
}
dump(RUN / "raw_model/pages_0001_0006_extraction.json", raw)

page_record_ids = {p: [] for p in range(1, 7)}
for rxn in reactions:
    for eid in rxn["evidence_ids"]:
        if ":page" in eid:
            p = int(eid.split(":page", 1)[1].split(":", 1)[0])
            if p in page_record_ids and rxn["reaction_id"] not in page_record_ids[p]:
                page_record_ids[p].append(rxn["reaction_id"])
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/abstract and graphical abstract", "finding": "Article metadata and a generic this-work reaction summary; no distinct concrete experiment.", "evidence_ids": E[1]},
            {"region": "two-column introduction", "finding": "Background and cited prior art discussion; no current-paper concrete experimental row.", "evidence_ids": E[1]}],
         "record_ids": page_record_ids[1], "unresolved": ["Graphical-abstract Ar/RF/tether variables and ligand L* remain generic."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Figure 1a", "finding": "Background bioactive structures; source_context, not current experiments.", "evidence_ids": E[2]},
            {"region": "Figure 1b", "finding": "Cited Cramer literature example, retained as literature_example.", "evidence_ids": E[2]},
            {"region": "Figure 1c", "finding": "This-work generic reaction scheme, not a concrete row.", "evidence_ids": E[2]},
            {"region": "Table 1 and footnotes", "finding": "Ten optimization rows; entry 9 duplicates Scheme 1 product 3a and is merged with page-3 evidence.", "evidence_ids": E[2]},
            {"region": "scope discussion prose", "finding": "Narrative cross-checks for 3a-u, including specific values for 3b, 3n, and 3o; no additional duplicate experiments.", "evidence_ids": E[2]}],
         "record_ids": page_record_ids[2], "unresolved": ["Table 1 footnote c does not explicitly restate scaled absolute charges/volumes; conflict with Scheme 1 is retained for entry 9."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Scheme 1 green general-reaction panel", "finding": "Generic CF3/CF2H imidoylation conditions and footnotes b-e.", "evidence_ids": E[3]},
            {"region": "Scheme 1 CF3 upper scope", "finding": "Products 3a-z plus scale demonstrations, each transcribed; 3a standard row deduplicated against Table 1 entry 9.", "evidence_ids": E[3]},
            {"region": "Scheme 1 CF3 isocyanide scope", "finding": "Products 3aa-ae transcribed as separate experiments.", "evidence_ids": E[3]},
            {"region": "Scheme 1 CF2H lower scope", "finding": "Products 4a-ab transcribed; 4a low-catalyst comparison is separate and 4y' is retained as a byproduct outcome.", "evidence_ids": E[3]},
            {"region": "X-ray renderings", "finding": "3u and 4a crystallographic renderings are characterization context, not new syntheses.", "evidence_ids": E[3]},
            {"region": "Scheme 1 footnotes", "finding": "Isolated-yield basis, catalyst-loading exceptions, 70 °C exception, and 0.1 mmol CF2H conditions transcribed.", "evidence_ids": E[3]}],
         "record_ids": page_record_ids[3], "unresolved": ["Most complete scope structures, exact absolute configurations, ee methods, and detailed workups require unavailable SI.", "Scale-up callouts do not restate all absolute charges; they were not calculated."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 2a", "finding": "Three allene-insertion cascade experiments 5a-c to 6a-c.", "evidence_ids": E[4]},
            {"region": "Scheme 2b", "finding": "Eight complex-molecule scope products 7a-h with yield and ee/de.", "evidence_ids": E[4]},
            {"region": "Scheme 2c", "finding": "Two distinct nitrations of 3a to 8a and 8b with acid volumes, room temperature, and 16/32 h.", "evidence_ids": E[4]},
            {"region": "X-ray renderings", "finding": "8a and 8b crystallographic renderings are characterization context.", "evidence_ids": E[4]},
            {"region": "body prose", "finding": "Cross-checks CF2H scope, 4y' byproduct, allene cascade, complex-molecule scope and nitration outcomes; no extra distinct experiment.", "evidence_ids": E[4]}],
         "record_ids": page_record_ids[4], "unresolved": ["Scheme 2b does not print absolute substrate charges; de major/minor configurations and selectivity methods are unreported in processed pages."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "Scheme 3a", "finding": "One true methyl-isocyanide control experiment yielding 10.", "evidence_ids": E[5]},
            {"region": "Scheme 3b", "finding": "Author-proposed favored/disfavored stereocontrol models; proposed_mechanism, not experimental proof.", "evidence_ids": E[5]},
            {"region": "summary prose", "finding": "Article-level summary only; no new experiment.", "evidence_ids": E[5]},
            {"region": "Supporting Information statement", "finding": "SI availability statement retained, but no SI file is linked or inferred.", "evidence_ids": E[5]},
            {"region": "author information/acknowledgments", "finding": "Administrative metadata; no chemical experiment.", "evidence_ids": E[5]}],
         "record_ids": page_record_ids[5], "unresolved": ["Control-product yield measurement basis and absolute configuration are not reported in processed pages."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "references 1-10", "finding": "Bibliographic citations only; no experiment performed in the current paper.", "evidence_ids": E[6]}],
         "record_ids": page_record_ids[6], "unresolved": []},
    ]
}
dump(RUN / "coverage/pages_0001_0006.json", coverage)

# Mark only this session's six page tasks as attempted; pages 7-8 remain pending.
tasks_path = RUN / "tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    if task["document_id"] == DOC and any(task["evidence_id"] == f"{DOC}:page{p}:text" for p in range(1, 7)):
        task["state"] = "attempted"
        task["reason"] = "Original page image inspected and mapped in coverage/pages_0001_0006.json; records remain needs_review."
dump(tasks_path, tasks)

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in reactions)}, ensure_ascii=False))
