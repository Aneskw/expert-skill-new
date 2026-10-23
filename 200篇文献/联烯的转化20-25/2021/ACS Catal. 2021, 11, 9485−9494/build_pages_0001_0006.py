#!/usr/bin/env python3
"""Build source-derived v3.1 records for PDF pages 1-6 only."""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st  # noqa: E402
from reaction_tools import dataset_audit  # noqa: E402


DOC = "doc_7aafc6041d0daf10"
TITLE = "Scope and Mechanistic Investigations of Pd-Catalyzed Coupling/Cyclization and Cycloisomerization of Allenyl Malonates"
RUN = ROOT / "runs/current"


def ev(page: int, kind: str = "image") -> str:
    return f"{DOC}:page{page}:{kind}"


def sf(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
       status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status)


def inferred(value, evidence, derivation):
    return sf(value, value, evidence, method="inference", derivation=derivation)


def set_fact(container, key, raw, normalized, evidence, method="original_extraction",
             derivation=None, status=None):
    container[key] = sf(raw, normalized, evidence, method, derivation, status)


class Builder:
    def __init__(self, bundle):
        self.bundle = bundle
        self.molecules = {}
        self.reactions = []

    def molecule(self, mid, label, name=None, evidence=(ev(2),), *, identity="source_identified",
                 chirality="unresolved", geometry="unresolved", representation="unresolved",
                 note=None, smiles=None, smiles_raw=None, axial=None, assignment=None,
                 mixture=None, double_bonds=None):
        if mid in self.molecules:
            m = self.molecules[mid]
            for eid in evidence:
                if not any(o["document_id"] == DOC and eid in o["evidence_ids"] for o in m["occurrences"]):
                    m["occurrences"].append({
                        "document_id": DOC,
                        "compound_label": sf(label, label, [eid]),
                        "evidence_ids": [eid],
                    })
            return m
        m = st.new_molecule(mid, DOC, list(evidence))
        m["identity_status"] = identity
        m["occurrences"] = [{
            "document_id": DOC,
            "compound_label": sf(label, label, evidence),
            "evidence_ids": list(evidence),
        }]
        m["names"]["original"] = sf(name or label, name or label, evidence)
        m["names"]["english"] = sf(name or label, name or label, evidence)
        m["names"]["abbreviations"] = sf()
        m["names"]["aliases"] = sf()
        m["cas_numbers"] = sf()
        if smiles:
            m["representation_type"] = "smiles"
            m["smiles"] = sf(smiles_raw or "connectivity interpreted from the source drawing", smiles,
                              evidence, method="structure_derivation",
                              derivation="Candidate connectivity was transcribed from the drawn structure; parsing does not verify source identity or allene axial stereochemistry.")
        else:
            m["representation_type"] = "unresolved"
            m["smiles"] = sf()
        stereo_basis = "Stereo state retained conservatively from the visible label/drawing; absence of a wedge or descriptor is not treated as achiral or racemic."
        m["stereochemistry"]["chirality_status"] = inferred(chirality, evidence, stereo_basis)
        m["stereochemistry"]["geometry_status"] = inferred(geometry, evidence, stereo_basis)
        m["stereochemistry"]["representation_status"] = inferred(representation, evidence, stereo_basis)
        m["stereochemistry"]["axial"] = sf(axial, axial, evidence) if axial is not None else sf()
        m["stereochemistry"]["assignment"] = sf(assignment, assignment, evidence) if assignment is not None else sf()
        m["stereochemistry"]["mixture"] = sf(mixture, mixture, evidence) if mixture is not None else sf()
        m["stereochemistry"]["double_bonds"] = sf(double_bonds, double_bonds, evidence) if double_bonds is not None else sf()
        m["stereochemistry"]["atom_centers"] = sf()
        m["stereochemistry"]["isomeric_smiles"] = sf()
        m["stereochemistry"]["identity_key"] = sf()
        m["identity_note"] = sf(note, note, evidence) if note else sf()
        m["external_properties"] = []
        self.molecules[mid] = m
        return m

    def participant(self, pid, mid, role, evidence, *, equiv=None, molpct=None, stage="reaction",
                    product_role=None):
        p = st.new_participant(pid, mid, role, list(evidence), stage=stage)
        if product_role:
            p["product_role"] = sf(product_role, product_role, evidence)
        if equiv is not None:
            p["charges"]["equivalents"] = sf(f"{equiv:g} equiv", {"number": equiv, "unit": "equiv"}, evidence)
            p["equivalence_reference"] = sf("relative to allenyl substrate", "m_2a_racemic", evidence)
        if molpct is not None:
            p["charges"]["mol_percent"] = sf(f"{molpct:g} mol%", {"number": molpct, "unit": "mol%"}, evidence)
            p["equivalence_reference"] = sf("relative to allenyl substrate", "m_2a_racemic", evidence)
        return p

    def base_reaction(self, rid, evidence, record_type, categories, location, name, transformation):
        r = st.new_reaction(f"{DOC}:{rid}", DOC, list(evidence))
        r["record_type"] = record_type
        r["experiment_categories"] = sf(categories, categories, evidence)
        r["location"] = sf(location, location, evidence)
        r["reaction_name"] = sf(name, name, evidence)
        r["transformation_type"] = inferred(transformation, evidence,
            "Transformation description is a normalized interpretation of the article heading, scheme caption, and visible reactant/product drawings.")
        r["review"] = {
            "status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
            "unresolved": ["No explicitly associated Supporting Information was supplied; unreported procedural and characterization details remain missing."],
        }
        r["completeness"] = "partial"
        self.reactions.append(r)
        return r

    def classify_none(self, r, evidence, reason):
        c = r["classification"]
        c["mapped_reaction_smiles"] = sf()
        c["bond_changes"] = []
        c["labels"] = sf([], [], evidence, method="inference", derivation=reason)
        c["coverage"] = sf("not applicable", "not_applicable", evidence, method="inference", derivation=reason)
        c["limitations"] = [reason]
        c["model_votes"] = []

    def classify_changes(self, r, evidence, changes, limitations):
        c = r["classification"]
        c["mapped_reaction_smiles"] = sf()
        facts = []
        expected = {}
        for ch in changes:
            deriv = ch.pop("derivation")
            facts.append(sf(ch.copy(), ch.copy(), evidence, method="structure_derivation", derivation=deriv))
            if ch["change_type"] in ("formed", "order_increased", "hydrogen_gain"):
                label = st.bond_label(ch["atoms"][0]["element"], ch["atoms"][1]["element"], ch["after_order"])
                expected.setdefault(label["code"], {**label, "supporting_change_ids": []})["supporting_change_ids"].append(ch["change_id"])
        c["bond_changes"] = facts
        c["labels"] = sf("source drawings interpreted under bond-formation-20261016", list(expected.values()), evidence,
                         method="structure_derivation",
                         derivation="Labels were derived from the explicitly listed partial net bond changes; no complete atom mapping was asserted.")
        c["coverage"] = sf("partial", "partial", evidence, method="inference",
                           derivation="Structures are generic, unresolved, or not fully atom-mapped; only unambiguous bond classes are retained.")
        c["limitations"] = limitations
        c["model_votes"] = []

    def add_step(self, r, evidence, operation, temperature=None, time=None, participants=None):
        s = st.new_step(f"step_{len(r['steps'])+1}", list(evidence))
        s["participant_ids"] = participants or [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
        s["operation"] = sf(operation, operation, evidence)
        if temperature is not None:
            raw, val = temperature
            s["temperature"] = sf(raw, val, evidence)
        if time is not None:
            raw, val = time
            s["time"] = sf(raw, val, evidence)
        r["steps"].append(s)
        return s

    def outcome(self, r, oid, kind, targets, evidence, raw, normalized, basis=None,
                product_role=None, yield_type=None, scope=None):
        o = st.new_outcome(oid, kind, targets, list(evidence))
        o["value"] = sf(raw, normalized, evidence)
        if basis is not None:
            o["measurement_basis"] = sf(basis, basis, evidence)
        if product_role is not None:
            o["product_role"] = sf(product_role, product_role, evidence)
        if yield_type is not None:
            o["yield_type"] = sf(yield_type, yield_type, evidence)
        if scope is not None:
            o["scope"] = sf(scope, scope, evidence)
        r["outcomes"].append(o)
        return o


def c_c_changes(prefix, count, target):
    out = []
    for i in range(count):
        out.append({
            "change_id": f"{prefix}_cc_{i+1}",
            "atoms": [{"atom_map": None, "element": "C", "site": f"{target}:siteA{i+1}"},
                      {"atom_map": None, "element": "C", "site": f"{target}:siteB{i+1}"}],
            "before_order": 0, "after_order": 1, "change_type": "formed", "target": target,
            "derivation": f"The {target} drawing shows a new carbon-carbon connection in the net cyclization/coupling; exact atom maps are unresolved.",
        })
    return out


def cycloisomerization_changes(prefix, target):
    return c_c_changes(prefix, 1, target) + [
        {"change_id": f"{prefix}_h_loss", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}],
         "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "target": target,
         "derivation": "A net C-H redistribution is required by the drawn allene-to-cyclopentene isomerization, but the donor/acceptor atom mapping is unresolved."},
        {"change_id": f"{prefix}_h_gain", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}],
         "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "target": target,
         "derivation": "A net C-H redistribution is required by the drawn allene-to-cyclopentene isomerization, but the donor/acceptor atom mapping is unresolved."},
    ]


def common_molecules(b: Builder):
    # Source-explicit names are kept, but CAS and unshown structures are not supplied from external knowledge.
    named = {
        "m_ph_i": ("PhI", "iodobenzene"), "m_pd2dba3": ("Pd2dba3", "tris(dibenzylideneacetone)dipalladium"),
        "m_dppm": ("dppm", "1,1-bis(diphenylphosphino)methane"), "m_dppe": ("dppe", "1,2-bis(diphenylphosphino)ethane"),
        "m_dppp": ("dppp", "1,3-bis(diphenylphosphino)propane"), "m_dppf": ("dppf", "1,1'-ferrocenediyl-bis(diphenylphosphine)"),
        "m_poph3": ("P(OPh)3", "triphenyl phosphite"), "m_pcf3ph3": ("P(p-CF3-Ph)3", "tris(4-trifluoromethylphenyl)phosphine"),
        "m_me_duphos": ("(R,R)-Me-DuPhos", "(R,R)-Me-DuPhos"), "m_tangphos": ("(S,S',R,R')-TangPhos", "TangPhos"),
        "m_josiphos": ("JOSIPHOS SL-J003-1", "JOSIPHOS SL-J003-1"), "m_walphos": ("WALPHOS SL-W003-1", "WALPHOS SL-W003-1"),
        "m_spanphos": ("(±)-SPANPHOS", "racemic SPANPHOS"), "m_lihmds": ("LiHMDS", "lithium bis(trimethylsilyl)amide"),
        "m_nahmds": ("NaHMDS", "sodium bis(trimethylsilyl)amide"), "m_khmds": ("KHMDS", "potassium bis(trimethylsilyl)amide"),
        "m_nah": ("NaH", "sodium hydride"), "m_kh": ("KH", "potassium hydride"), "m_k2co3": ("K2CO3", "potassium carbonate"),
        "m_tbab": ("TBAB", "tetrabutylammonium bromide"), "m_pdpph34": ("Pd(PPh3)4", "tetrakis(triphenylphosphine)palladium(0)"),
        "m_dmf": ("DMF", "dimethylformamide"), "m_thf": ("THF", "tetrahydrofuran"), "m_dmso": ("DMSO", "dimethyl sulfoxide"),
        "m_mecn": ("MeCN", "acetonitrile"), "m_rucl3": ("RuCl3", "ruthenium(III) chloride"),
        "m_naio4": ("NaIO4", "sodium periodate"), "m_cecl3_7h2o": ("CeCl3·7H2O", "cerium(III) chloride heptahydrate"),
        "m_meoh": ("MeOH", "methanol"), "m_h2o": ("H2O", "water"), "m_etoac": ("EtOAc", "ethyl acetate"),
        "m_etona": ("EtONa", "sodium ethoxide"), "m_etoh": ("EtOH", "ethanol"), "m_pdpph34_lit": ("Pd(PPh3)4", "tetrakis(triphenylphosphine)palladium(0)"),
    }
    for mid, (label, name) in named.items():
        chirality = "racemic" if mid == "m_spanphos" else "unresolved"
        b.molecule(mid, label, name, [ev(3), ev(4), ev(5)] if mid in {"m_pd2dba3", "m_dppm", "m_dppf", "m_nahmds", "m_kh", "m_dmf", "m_ph_i"} else [ev(2), ev(3), ev(5)], chirality=chirality,
                   note="Named reaction component; a complete source-verified structure was not encoded from these bounded pages.")
    smiles_2a = "CCCCCC=C=CCC(C(=O)OC)C(=O)OC"
    b.molecule("m_2a_racemic", "(±)-2a", "racemic allenyl malonate 2a", [ev(3)],
               chirality="racemic", geometry="unresolved", representation="partial", smiles=smiles_2a,
               smiles_raw="drawn connectivity of (±)-2a in Table 1",
               axial={"axis": "allene", "assignment": "racemic; individual enantiomer configuration not encoded"},
               mixture={"type": "racemate", "ratio": "not numerically stated"},
               note="The allene is explicitly marked (±); RDKit SMILES does not encode the source allene axial configuration.")
    b.molecule("m_2a_Ra_86ee", "(R_a)-2a", "(R_a)-allenyl malonate 2a, 86% ee", [ev(5)],
               chirality="specified", geometry="unresolved", representation="not_encoded", smiles=smiles_2a,
               smiles_raw="drawn connectivity of (R_a)-2a in Scheme 3",
               axial={"axis": "allene", "configuration": "R_a", "ee": "86%"},
               assignment={"absolute_configuration": "R_a", "enantiomeric_excess": "86%", "method": "not reported on pages 1-6"},
               note="Allene axial configuration is source-explicit but not encoded in the candidate SMILES.")
    for mid, label, note, chiral, geom in [
        ("m_3a", "3a", "cyclopentene product; stereochemistry is not fully assigned in Table 1/Table 2", "unspecified", "unresolved"),
        ("m_4a_table", "4a", "vinyl cyclopropane product in Table 1; E/Z and absolute configuration are not stated there", "unspecified", "unspecified"),
        ("m_5", "5", "cycloisomerized cyclopentene product", "unspecified", "unresolved"),
        ("m_6", "6", "dihydroxylated cyclopentane product reported as a 1.4:1 diastereomeric mixture", "mixture", "not_applicable"),
    ]:
        b.molecule(mid, label, label, [ev(3), ev(4)] if mid in {"m_3a", "m_5"} else [ev(3)] if mid == "m_4a_table" else [ev(4)],
                   chirality=chiral, geometry=geom, representation="unresolved", note=note,
                   mixture={"type": "diastereomeric mixture", "ratio": "1.4:1"} if mid == "m_6" else None)


def add_literature_and_generic(b: Builder):
    e = [ev(2)]
    generic_mols = [
        ("m_lit_cazes_allene", "Cazes allene", "generic allene with C15H7 substituent"),
        ("m_lit_cazes_malonate", "Cazes malonate", "generic malonate E1/E2/R"),
        ("m_lit_cazes_products", "Cazes product mixture", "generic olefin/dendralene mixture; poor yields and E/Z ratios"),
        ("m_lit_ma_substrate", "Ma allenyl malonate", "generic allenyl malonate A/B"),
        ("m_lit_ma_product_A", "Ma cyclopentene product", "favored product under pathway A conditions"),
        ("m_lit_ma_product_B", "Ma vinyl cyclopropane product", "favored product under pathway B conditions"),
        ("m_generic_allenyl_malonate_1", "1", "generic allenyl lactone substrate with R/E variables"),
        ("m_generic_allenyl_malonate_2", "2", "generic acyclic allenyl malonate substrate with R/E variables"),
        ("m_generic_allenyl_substrate_family", "1 or 2", "generic alternative allenyl malonate/lactone substrate family"),
        ("m_generic_cross_products", "generic 3/4 products", "generic cyclopentene or cyclopropane products"),
        ("m_generic_cycloisomer_product", "generic cycloisomerization product", "generic cyclopentene product without Ar-X"),
    ]
    for mid, label, note in generic_mols:
        b.molecule(mid, label, label, e, identity="unresolved", chirality="unresolved", geometry="unresolved", representation="unresolved", note=note)
    records = []
    r = b.base_reaction("scheme1_cazes_literature", e, "literature_example", ["other"], "Scheme 1, Cazes arrow",
                        "Cazes intermolecular malonate/allene reaction", "literature-reported intermolecular coupling")
    r["participants"] = [b.participant("allene", "m_lit_cazes_allene", "reactant", e), b.participant("malonate", "m_lit_cazes_malonate", "reactant", e),
                         b.participant("pd", "m_pd2dba3", "catalyst", e), b.participant("dppe", "m_dppe", "ligand", e),
                         b.participant("etona", "m_etona", "reagent", e), b.participant("etoh", "m_etoh", "solvent", e),
                         b.participant("products", "m_lit_cazes_products", "product", e, product_role="reported mixture")]
    b.add_step(r, e, "literature-example reaction arrow", participants=[p["participant_id"] for p in r["participants"]])
    b.classify_changes(r, e, c_c_changes("cazes", 1, "generic product mixture"), ["Prior work cited in Scheme 1; generic R/E variables and multiple products preclude full mapping."])
    records.append(r["reaction_id"])
    for suffix, product, conditions in [
        ("A", "m_lit_ma_product_A", "Pd(PPh3)4; 1.1 equiv NaH; DMSO; 80 °C"),
        ("B", "m_lit_ma_product_B", "Pd(PPh3)4; K2CO3; CH3CN; reflux"),
    ]:
        r = b.base_reaction(f"scheme1_ma_{suffix.lower()}_literature", e, "literature_example", ["other"], f"Scheme 1, Ma arrow {suffix}",
                            "Ma Pd-catalyzed coupling/cyclization", "literature-reported regioselective coupling/cyclization")
        r["participants"] = [b.participant("substrate", "m_lit_ma_substrate", "reactant", e), b.participant("phi", "m_ph_i", "reactant", e),
                             b.participant("catalyst", "m_pdpph34_lit", "catalyst", e), b.participant("product", product, "product", e, product_role="favored product")]
        if suffix == "A":
            r["participants"] += [b.participant("nah", "m_nah", "reagent", e), b.participant("solvent", "m_dmso", "solvent", e)]
            b.add_step(r, e, conditions, ("80 °C", {"number": 80, "unit": "°C"}))
        else:
            r["participants"] += [b.participant("k2co3", "m_k2co3", "reagent", e), b.participant("solvent", "m_mecn", "solvent", e)]
            b.add_step(r, e, conditions, ("reflux", "reflux"))
        b.classify_changes(r, e, c_c_changes(f"ma_{suffix}", 2, f"Ma product {suffix}"), ["Prior work cited in Scheme 1; generic substituents preclude complete mapping."])
        records.append(r["reaction_id"])
    r = b.base_reaction("scheme1_current_cross_coupling_generic", e, "generic_scheme", ["other"], "Scheme 1, This work, Ar-X arrow",
                        "Pd-catalyzed cross-coupling/cyclization of allenyl malonates", "cross-coupling/cyclization to cyclopentenes or cyclopropanes")
    r["participants"] = [b.participant("substrate_family", "m_generic_allenyl_substrate_family", "reactant", e),
                         b.participant("arx", "m_ph_i", "reactant", e), b.participant("pd", "m_pd2dba3", "catalyst", e, molpct=5),
                         b.participant("dppm", "m_dppm", "ligand", e, molpct=10), b.participant("nahmds", "m_nahmds", "reagent", e, equiv=1.0),
                         b.participant("dmf", "m_dmf", "solvent", e), b.participant("products", "m_generic_cross_products", "product", e)]
    for p in r["participants"]:
        if p["charges"]["equivalents"]["normalized_value"] is not None or p["charges"]["mol_percent"]["normalized_value"] is not None:
            p["equivalence_reference"] = sf("relative to the chosen allenyl substrate 1 or 2", "m_generic_allenyl_substrate_family", e)
    r["review"]["unresolved"].append("Scheme 1 is a generic summary: Ar-X, E, and R groups are variable and the 1:2 Pd2dba3:dppm notation is not a complete procedure.")
    b.add_step(r, e, "generic reaction arrow", ("rt", "rt"))
    b.classify_changes(r, e, c_c_changes("current_generic", 2, "generic cross-coupled cyclic product"), ["Generic scheme; exact atom mapping and product branch are unresolved."])
    records.append(r["reaction_id"])
    r = b.base_reaction("scheme1_current_cycloisomerization_generic", e, "generic_scheme", ["other"], "Scheme 1, This work, no Ar-X arrow",
                        "Pd-catalyzed cycloisomerization of allenyl malonates", "cycloisomerization to cyclopentene")
    r["participants"] = [b.participant("sub", "m_generic_allenyl_malonate_2", "reactant", e), b.participant("pd", "m_pd2dba3", "catalyst", e),
                         b.participant("poph3", "m_poph3", "ligand", e), b.participant("kh", "m_kh", "reagent", e),
                         b.participant("dmf", "m_dmf", "solvent", e), b.participant("product", "m_generic_cycloisomer_product", "product", e)]
    b.add_step(r, e, "generic reaction arrow without Ar-X", ("rt", "rt"))
    b.classify_changes(r, e, cycloisomerization_changes("cycloisom_generic", "generic cyclopentene"), ["Generic scheme; carbon and hydrogen atom mapping is partial."])
    records.append(r["reaction_id"])
    return records


def add_table1(b: Builder):
    rows = [
        (1,"m_dppm","m_dmf","m_lihmds",1.2,"rt",{"3a":10,"2a":49}),
        (2,"m_dppm","m_dmf","m_nahmds",1.2,"rt",{"3a":65}),
        (3,"m_dppm","m_dmf","m_nahmds",1.0,"rt",{"3a":71}),
        (4,"m_dppm","m_dmf","m_khmds",1.2,"rt",{"3a":33}),
        (5,"m_dppm","m_dmf","m_nah",1.2,"rt",{"3a":71,"3a_iso":67,"5":0}),
        (6,"m_dppe","m_dmf","m_nah",1.2,"rt",{"3a":38,"5":7}),
        (7,"m_dppp","m_dmf","m_nah",1.2,"rt",{"3a":51,"5":6}),
        (8,"m_dppf","m_dmf","m_nah",1.2,"rt",{"3a":44}),
        (9,"m_dppf","m_dmf","m_nah",1.2,"66 °C",{"3a":54}),
        (10,"m_dppf","m_thf","m_nah",1.2,"rt",{"3a":0,"2a":85}),
        (11,"m_dppf","m_thf","m_nah",1.2,"reflux",{"3a":52}),
        (12,"m_dppf","m_dmso","m_nah",1.2,"rt",{"3a":"trace","2a":70}),
        (13,"m_dppf","m_dmf","m_kh",1.2,"rt",{"3a":0,"2a":44,"4a":34}),
        (14,"m_poph3","m_dmf","m_nah",1.2,"rt",{"3a":16,"2a":0,"4a":0,"5":24}),
        (15,"m_pcf3ph3","m_dmf","m_nah",1.2,"rt",{"3a":25,"2a":0,"4a":0,"5":12}),
        (16,"m_me_duphos","m_dmf","m_nah",1.2,"rt",{"3a":6,"2a":83,"4a":0,"5":0}),
        (17,"m_tangphos","m_dmf","m_nah",1.2,"rt",{"3a":0,"2a":93,"4a":0,"5":0}),
        (18,"m_josiphos","m_dmf","m_nah",1.2,"rt",{"3a":0,"2a":83,"4a":0,"5":0}),
        (19,"m_walphos","m_dmf","m_nah",1.2,"rt",{"3a":32,"2a":18,"4a":0,"5":0}),
        (20,"m_spanphos","m_dmf","m_nah",1.2,"rt",{"3a":58,"2a":8,"4a":0,"5":0}),
    ]
    ids=[]
    for n, ligand, solvent, base, base_eq, temp, results in rows:
        e=[ev(3)]
        r=b.base_reaction(f"table1_entry_{n:02d}",e,"specific_experiment",["optimization"],f"Table 1, entry {n}",
                          "Selected reaction optimization for cross-coupling/cyclization of 2a to 3a","Pd-catalyzed coupling/cyclization optimization")
        parts=[b.participant("substrate","m_2a_racemic","reactant",e),b.participant("phi","m_ph_i","reactant",e,equiv=1.2),
               b.participant("pd","m_pd2dba3","catalyst",e,molpct=5),b.participant("ligand",ligand,"ligand",e,molpct=10),
               b.participant("base",base,"reagent",e,equiv=base_eq),b.participant("solvent",solvent,"solvent",e)]
        if "3a" in results or "3a_iso" in results: parts.append(b.participant("product_3a","m_3a","product",e,product_role="target cyclopentene"))
        if "4a" in results: parts.append(b.participant("product_4a","m_4a_table","product",e,product_role="vinyl cyclopropane"))
        if "5" in results: parts.append(b.participant("product_5","m_5","product",e,product_role="cycloisomerization product"))
        r["participants"]=parts
        r["solvent_conditions"]["composition"]=sf(b.molecules[solvent]["names"]["original"]["raw_value"],b.molecules[solvent]["names"]["original"]["normalized_value"],e)
        tval={"number":66,"unit":"°C"} if temp=="66 °C" else temp
        b.add_step(r,e,"Table 1 optimization reaction",(temp,tval))
        r["procedure"]["original_text"]=sf({"header":"5 mol% Pd2dba3; 10 mol% ligand; 1.2 equiv PhI; base; solvent; temperature","entry":n},
                                                   f"Table 1 entry {n}",e)
        changes=[]
        if isinstance(results.get("3a"),int) and results["3a"]>0 or results.get("3a")=="trace": changes+=c_c_changes(f"t1e{n}_3a",2,"3a")
        if isinstance(results.get("4a"),int) and results["4a"]>0: changes+=c_c_changes(f"t1e{n}_4a",2,"4a")
        if isinstance(results.get("5"),int) and results["5"]>0: changes+=cycloisomerization_changes(f"t1e{n}_5","5")
        for key,val in results.items():
            if key=="3a_iso":
                b.outcome(r,"yield_3a_isolated","yield",["m_3a"],e,"67% (isolated)",{"number":67,"unit":"%"},"isolated yield",product_role="target product",yield_type="isolated")
            elif key=="2a":
                b.outcome(r,"recovery_2a","recovery",["m_2a_racemic"],e,f"{val}%",{"number":val,"unit":"%"},"1H NMR with Ph3SiMe as internal standard",product_role="recovered starting material")
            else:
                mid={"3a":"m_3a","4a":"m_4a_table","5":"m_5"}[key]
                if val=="trace":
                    b.outcome(r,f"trace_{key}","trace",[mid],e,"trace","trace","1H NMR with Ph3SiMe as internal standard",product_role="observed product")
                else:
                    b.outcome(r,f"yield_{key}","yield",[mid],e,f"{val}%",{"number":val,"unit":"%"},"1H NMR with Ph3SiMe as internal standard",product_role="reported product",yield_type="NMR")
        if changes:
            b.classify_changes(r,e,changes,["Partial product-specific net changes were read from the drawings; exact atom mapping is not asserted.","Rows with multiple products combine product-specific changes in one experiment record."])
        else:
            b.classify_none(r,e,"No positive product formation is reported in this optimization row; zero yields do not justify an intended-product bond classification.")
        r["review"]["unresolved"] += ["Reaction time, scale, concentration, atmosphere, pressure, workup, and purification are not reported on page 3.","Blank Table 1 cells were preserved as unreported rather than converted to zero."]
        ids.append(r["reaction_id"])
    return ids


def add_scope_and_scaleup(b: Builder):
    data=[("3a",67,None,None),("3b",73,None,None),("3c",50,None,None),("3d",72,"average NMR yield over two runs",None),
          ("3e",67,None,None),("3f",73,None,None),("3g",77,None,None),("3h",51,None,None),("3i",60,None,None),
          ("3j",70,None,None),("3k",65,None,None),("3l",39,None,None),("3m",66,None,None),("3n",87,None,None),
          ("3o",48,None,None),("3p",38,"yield based on isolated mixture of 3p/2a",None),("3q",52,"NMR yield","1.2:1 dr"),
          ("3r",60,None,"1.3:1 dr"),("3s",59,None,"1:1 dr")]
    aryl={"3a":"iodobenzene","3b":"4-iodobenzotrifluoride","3c":"4-iodobenzaldehyde","3d":"4-iodoacetophenone",
          "3e":"ethyl 4-iodobenzoate","3f":"1-bromo-4-iodobenzene","3g":"3-iodopyridine","3h":"1-iodonaphthalene",
          "3i":"1-iodo-3,5-dimethylbenzene","3j":"1-iodo-2-methylbenzene","3k":"1-iodo-2-methoxybenzene",
          "3l":"chloro-nitro iodoarene (positions unresolved from bounded text)","3m":"2-iodothiophene"}
    ids=[]
    for label,yield_n,basis,dr in data:
        e=[ev(4)]
        pmid=f"m_{label}"
        if pmid not in b.molecules:
            b.molecule(pmid,label,label,e,chirality="mixture" if dr else "unspecified",geometry="unresolved",representation="unresolved",
                       mixture={"type":"diastereomeric mixture","ratio":dr.replace(" dr","")} if dr else None,
                       note=f"Product {label} as drawn in Table 2; complete source-verified SMILES and absolute stereochemistry are not available from the bounded pages.")
        if label in [f"3{x}" for x in "abcdefghijklm"]:
            sub="m_2a_racemic"
        else:
            sub=f"m_precursor_{label}"
            precursor_label="2"+label[1:]
            b.molecule(sub,precursor_label,f"allenyl malonate/lactone {precursor_label}, precursor to {label}",e,chirality="unspecified",geometry="unresolved",representation="unresolved",
                       note="The scope drawing is legible as a specific connectivity but is retained unresolved rather than assigning an unchecked SMILES.")
        amid=f"m_aryl_{label}"
        aname=aryl.get(label,"iodobenzene")
        b.molecule(amid,"ArI",aname,e,chirality="unresolved",geometry="unresolved",representation="unresolved",note=f"Coupling partner for {label}: {aname}.")
        r=b.base_reaction(f"table2_{label}",e,"specific_experiment",["substrate_scope"],f"Table 2, product {label}",
                          "Scope of allenyl malonate cross-coupling/cyclization","Pd-catalyzed cross-coupling/cyclization")
        r["general_scheme_ids"]=[f"{DOC}:scheme1_current_cross_coupling_generic"]
        r["participants"]=[b.participant("substrate",sub,"reactant",e),b.participant("aryl_iodide",amid,"reactant",e,equiv=1.2),
                           b.participant("pd","m_pd2dba3","catalyst",e,molpct=5),b.participant("dppm","m_dppm","ligand",e,molpct=10),
                           b.participant("nahmds","m_nahmds","reagent",e,equiv=1.0),b.participant("dmf","m_dmf","solvent",e),
                           b.participant("product",pmid,"product",e,product_role="scope product")]
        # Fix denominator for scope entries whose substrate is not 2a.
        for p in r["participants"]:
            if p["charges"]["equivalents"]["normalized_value"] is not None or p["charges"]["mol_percent"]["normalized_value"] is not None:
                p["equivalence_reference"]=sf("relative to allenyl substrate",sub,e)
        r["solvent_conditions"]["composition"]=sf("DMF","DMF",e)
        temp="66 °C" if label in {"3p","3s"} else "rt"
        tval={"number":66,"unit":"°C"} if temp=="66 °C" else "rt"
        b.add_step(r,e,"Table 2 scope reaction",(temp,tval),("overnight","overnight"))
        if label=="3p": r["special_conditions"]["reactor"]=sf("pressure tube","pressure tube",e)
        b.outcome(r,"yield","yield",[pmid],e,f"{yield_n}%",{"number":yield_n,"unit":"%"},basis,product_role="scope product",
                  yield_type="NMR" if basis and "NMR" in basis else None,scope=f"allenyl malonate cross-coupling/cyclization to {label}")
        if dr:
            ratio=dr.replace(" dr","")
            b.outcome(r,"dr","selectivity",[pmid],e,dr,{"metric":"dr","ratio":ratio,"isomer_order":"major:minor; configurations unassigned"},None,
                      product_role="diastereomeric product mixture",scope="diastereomeric ratio; measurement method not reported on page 4")
        b.classify_changes(r,e,c_c_changes(f"scope_{label}",2,label),["Two net C-C formations are evident for coupling plus ring closure, but exact atom mapping and hydrogen changes are unresolved."])
        r["review"]["unresolved"] += ["Unmarked Table 2 yield measurement basis is not stated on page 4; only footnoted NMR/mixture bases were populated.","Reaction scale, concentration, atmosphere, workup, purification, and characterization are not reported on page 4."]
        ids.append(r["reaction_id"])
    # Scale-up to 3r.
    e=[ev(4)]
    r=b.base_reaction("scheme_table2_scaleup_3r",e,"specific_experiment",["scale_up"],"Table 2, Scale up and postfunctionalization, first arrow",
                      "2 mmol scale synthesis of 3r","Pd-catalyzed cross-coupling/cyclization scale-up")
    sub="m_precursor_3r"
    r["participants"]=[b.participant("substrate",sub,"reactant",e),b.participant("phi","m_aryl_3r","reactant",e,equiv=1.2),
                       b.participant("pd","m_pd2dba3","catalyst",e,molpct=5),b.participant("dppm","m_dppm","ligand",e,molpct=10),
                       b.participant("nahmds","m_nahmds","reagent",e,equiv=1.0),b.participant("dmf","m_dmf","solvent",e),
                       b.participant("product","m_3r","product",e,product_role="scale-up product")]
    for p in r["participants"]:
        if p["charges"]["equivalents"]["normalized_value"] is not None or p["charges"]["mol_percent"]["normalized_value"] is not None:
            p["equivalence_reference"]=sf("relative to allenyl substrate",sub,e)
    r["reaction_scale"]=sf("2 mmol scale",{"reference_molecule_id":sub,"amount":{"number":2,"unit":"mmol"}},e)
    r["solvent_conditions"]["composition"]=sf("DMF","DMF",e)
    b.add_step(r,e,"2 mmol scale cross-coupling/cyclization",("rt","rt"),("overnight","overnight"))
    b.outcome(r,"yield","yield",["m_3r"],e,"60%",{"number":60,"unit":"%"},None,product_role="scale-up product")
    b.outcome(r,"dr","selectivity",["m_3r"],e,"dr 1.3:1",{"metric":"dr","ratio":"1.3:1","isomer_order":"major:minor; configurations unassigned"},None,product_role="diastereomeric mixture")
    b.classify_changes(r,e,c_c_changes("scaleup_3r",2,"3r"),["Partial net C-C changes; no complete atom mapping."])
    r["review"]["unresolved"] += ["Yield measurement basis, concentration, atmosphere, pressure, workup, and purification are not reported on page 4."]
    ids.append(r["reaction_id"])
    # Dihydroxylation to 6.
    r=b.base_reaction("scheme_table2_dihydroxylation_6",e,"specific_experiment",["derivatization"],"Table 2, Scale up and postfunctionalization, second arrow",
                      "RuCl3-catalyzed dihydroxylation of 3r","alkene dihydroxylation")
    r["participants"]=[b.participant("substrate","m_3r","reactant",e),b.participant("rucl3","m_rucl3","catalyst",e,molpct=3.75),
                       b.participant("naio4","m_naio4","reagent",e,equiv=2.5),b.participant("cecl3","m_cecl3_7h2o","additive",e,molpct=30),
                       b.participant("meoh","m_meoh","solvent",e),b.participant("h2o","m_h2o","solvent",e),b.participant("etoac","m_etoac","solvent",e),
                       b.participant("product","m_6","product",e,product_role="dihydroxylated cyclopentane")]
    for p in r["participants"]:
        if p["charges"]["equivalents"]["normalized_value"] is not None or p["charges"]["mol_percent"]["normalized_value"] is not None:
            p["equivalence_reference"]=sf("relative to 3r","m_3r",e)
    r["solvent_conditions"]["composition"]=sf("4:1:2 MeCN:H2O:EtOAc",[{"solvent":"MeCN","parts":4},{"solvent":"H2O","parts":1},{"solvent":"EtOAc","parts":2}],e)
    r["solvent_conditions"]["mixture_ratio"]=sf("4:1:2 MeCN:H2O:EtOAc","4:1:2",e)
    b.add_step(r,e,"initial low-temperature stage; separate duration not reported",("0 °C",{"number":0,"unit":"°C"}))
    b.add_step(r,e,"warm from 0 °C to room temperature; the arrow gives 2.5 h for the overall 0 °C-to-rt operation",("rt","rt"),("2.5 h",{"number":2.5,"unit":"h"}))
    b.outcome(r,"yield","yield",["m_6"],e,"68%",{"number":68,"unit":"%"},None,product_role="dihydroxylated product")
    b.outcome(r,"dr","selectivity",["m_6"],e,"1.4:1 dr",{"metric":"dr","ratio":"1.4:1","isomer_order":"major:minor; configurations unassigned"},None,product_role="diastereomeric product mixture")
    changes=[]
    for i in (1,2):
        changes.append({"change_id":f"dihydrox_c_o_{i}","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"O"}],"before_order":0,"after_order":1,"change_type":"formed","target":"6",
                        "derivation":"The product drawing shows installation of two C-O bonds across the alkene; exact oxygen atom sources and atom maps are unresolved."})
        changes.append({"change_id":f"dihydrox_o_h_{i}","atoms":[{"atom_map":None,"element":"O"},{"atom_map":None,"element":"H"}],"before_order":0,"after_order":1,"change_type":"hydrogen_gain","target":"6",
                        "derivation":"The product drawing contains two O-H bonds; exact hydrogen source and atom maps are unresolved."})
    changes.append({"change_id":"dihydrox_cc_decrease","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"C"}],"before_order":2,"after_order":1,"change_type":"order_decreased","target":"6",
                    "derivation":"The alkene C=C becomes a C-C single bond in the diol product; this decrease is retained but not labeled as a new C-C sigma connection."})
    b.classify_changes(r,e,changes,["C-O and O-H formations are source-structure derived; oxygen/hydrogen sources and complete atom mapping remain unresolved."])
    r["review"]["unresolved"] += ["Yield and dr measurement methods, addition order, atmosphere, workup, and purification are not reported on page 4.","The arrow reports a 0 °C-to-rt trajectory over 2.5 h but does not split separate addition, ramp, and hold durations."]
    ids.append(r["reaction_id"])
    return ids


def add_mechanism_and_chirality(b: Builder):
    ids=[]
    e=[ev(5)]
    for label in ["2.1","2.2","2.3a","2.3b","2.4a","2.4b","2.6","generic 3","generic 4"]:
        b.molecule("m_mech_"+label.replace(".","_").replace(" ","_"),label,label,e,identity="unresolved",chirality="unresolved",geometry="unresolved",representation="unresolved",
                   note="Generic mechanistic entity containing R/Ar/E and/or PdLn variables; not a concrete source-verified molecule.")
    r=b.base_reaction("scheme2_proposed_mechanism",e,"proposed_mechanism",["mechanistic_study"],"Scheme 2",
                      "Potential mechanisms for Pd-catalyzed cross-coupling/cyclization","author-proposed catalytic pathways A-D")
    r["participants"]=[b.participant(f"entity_{i}","m_mech_"+label.replace(".","_").replace(" ","_"),"unknown",e) for i,label in enumerate(["2.1","2.2","2.3a","2.3b","2.4a","2.4b","2.6","generic 3","generic 4"],1)]
    inventory=["oxidative addition of Ar-I to Pd(0)","coordination/deprotonation to 2.2","migratory insertion to 2.3a","sigma-pi-sigma interconversion through 2.3a/2.3b and 2.4a/2.4b","Path A to cyclopropane 4","Path B to cyclopentene 3 with racemization","Path C nucleopalladation to 2.6 then product 3 with chirality transfer","Path D direct cyclization from 2.3a"]
    r["notes"]=[sf(item,item,e) for item in inventory]
    b.classify_none(r,e,"Scheme 2 is an author-proposed mechanism with generic intermediates, not an observed net experiment; no net classification is asserted.")
    r["review"]["unresolved"]=["Generic variables and metal-bound intermediates are not encoded as unique SMILES.","Competing mechanistic pathways are hypotheses; the record does not promote them to experimentally verified elementary steps."]
    ids.append(r["reaction_id"])
    # Scheme 3 product entities kept distinct because the article says the new KH/dppf conditions give the opposite stereoisomer to Ma's.
    b.molecule("m_3a_chiral_84ee","3a (84% ee)","enantioenriched cyclopentene 3a",e,chirality="partially_specified",geometry="unresolved",representation="not_encoded",
               assignment={"ee":"84%","absolute_configuration":"not assigned on page 5"},note="Enantioenriched product; absolute central configuration is not stated on page 5.")
    for mid,label,geom,ee,note in [
        ("m_4a_ma_E","4a E isomer under Ma conditions","E","12%","Absolute configuration not assigned; ee reported for the E isomer."),
        ("m_4a_ma_Z","4a Z isomer under Ma conditions","Z","−3%","Absolute configuration not assigned; signed ee reported for the Z isomer."),
        ("m_4a_new_E","4a major E isomer under new conditions","E","87%","Article prose states the stereoisomer is opposite to that formed under Ma conditions; absolute label is not given."),
        ("m_4a_new_Z","4a minor Z isomer under new conditions","Z",None,"Minor Z isomer; ee not reported on page 5."),
    ]:
        b.molecule(mid,label,label,e,chirality="partially_specified" if ee else "unspecified",geometry="specified",representation="not_encoded",
                   assignment={"ee":ee,"absolute_configuration":"unassigned"} if ee else None,
                   double_bonds={"configuration":geom,"basis":"Scheme 3 E/Z label"},note=note)
    scheme3=[
        ("scheme3a_chirality_transfer_3a","Our new conditions",[("product","m_3a_chiral_84ee")],
         [("m_pd2dba3",5,"catalyst"),("m_dppm",10,"ligand")],[("m_ph_i",1.20,"reactant"),("m_nahmds",1.0,"reagent")],"m_dmf","rt",
         [("ee_3a","selectivity",["m_3a_chiral_84ee"],"84% ee",{"metric":"ee","number":84,"unit":"%","configuration":"unassigned"}),
          ("fidelity","selectivity",["m_3a_chiral_84ee"],"98% fidelity",{"metric":"chirality-transfer fidelity","number":98,"unit":"%"})]),
        ("scheme3b_ma_conditions_4a","Ma's conditions",[("product_E","m_4a_ma_E"),("product_Z","m_4a_ma_Z")],
         [("m_pdpph34",10,"catalyst"),("m_tbab",5,"additive")],[("m_ph_i",1.20,"reactant"),("m_k2co3",4.0,"reagent")],"m_mecn","reflux",
         [("ee_E","selectivity",["m_4a_ma_E"],"12% ee_E",{"metric":"ee","number":12,"unit":"%","isomer":"E","configuration":"unassigned"}),
          ("ee_Z","selectivity",["m_4a_ma_Z"],"−3% ee_Z",{"metric":"signed ee","signed_percent":-3,"isomer":"Z","configuration/order":"source sign retained; absolute assignment unreported"}),
          ("ez","selectivity",["m_4a_ma_E","m_4a_ma_Z"],"1:1.4 E:Z",{"metric":"E:Z","ratio":"1:1.4","order":"E:Z"})]),
        ("scheme3c_new_conditions_4a","Our new conditions",[("product_E","m_4a_new_E"),("product_Z","m_4a_new_Z")],
         [("m_pd2dba3",5,"catalyst"),("m_dppf",10,"ligand")],[("m_ph_i",1.20,"reactant"),("m_kh",1.0,"reagent")],"m_dmf","rt",
         [("ee_E","selectivity",["m_4a_new_E"],"87% ee_E",{"metric":"ee","number":87,"unit":"%","isomer":"E","configuration":"unassigned"}),
          ("ez","selectivity",["m_4a_new_E","m_4a_new_Z"],">20:1 E:Z",{"metric":"E:Z","ratio":">20:1","order":"E:Z","qualifier":"gt"})]),
    ]
    for rid,cond,products,loaded,equivs,solvent,temp,outs in scheme3:
        r=b.base_reaction(rid,e,"specific_experiment",["mechanistic_study"],f"Scheme 3, {rid[7]}","Transfer-of-chirality experiment","axial-to-center chirality-transfer coupling/cyclization")
        r["participants"]=[b.participant("substrate","m_2a_Ra_86ee","reactant",e)]
        for pid,mid in products: r["participants"].append(b.participant(pid,mid,"product",e,product_role="chirality-transfer product"))
        for mid,pct,role in loaded: r["participants"].append(b.participant(mid.replace("m_",""),mid,role,e,molpct=pct))
        for mid,eq,role in equivs: r["participants"].append(b.participant(mid.replace("m_",""),mid,role,e,equiv=eq))
        r["participants"].append(b.participant("solvent",solvent,"solvent",e))
        # Denominator is the enantioenriched allene in Scheme 3.
        for p in r["participants"]:
            if p["charges"]["equivalents"]["normalized_value"] is not None or p["charges"]["mol_percent"]["normalized_value"] is not None:
                p["equivalence_reference"]=sf("relative to (R_a)-2a","m_2a_Ra_86ee",e)
        r["solvent_conditions"]["composition"]=sf(b.molecules[solvent]["names"]["original"]["raw_value"],b.molecules[solvent]["names"]["original"]["normalized_value"],e)
        b.add_step(r,e,cond,(temp,temp),("overnight","overnight"))
        r["notes"].append(sf("(R_a)-2a, 86% ee","Input allene (R_a)-2a has 86% ee",e))
        for oid,kind,targets,raw,norm in outs:
            b.outcome(r,oid,kind,targets,e,raw,norm,None,product_role="stereochemical outcome",scope="metric and isomer order retained exactly; measurement method not reported on page 5")
        changes=c_c_changes(rid,2,"3a" if "3a" in rid else "4a")
        b.classify_changes(r,e,changes,["Two net C-C formations are evident from the coupling/cyclization drawings; no complete atom mapping or absolute product configuration is asserted."])
        r["review"]["unresolved"] += ["Scheme 3 reports no product yield, reaction scale, concentration, atmosphere, workup, purification, or selectivity measurement method.","The allene R_a configuration is source-explicit but not encoded in the candidate SMILES; product absolute configurations are unassigned on page 5."]
        ids.append(r["reaction_id"])
    return ids


def add_page6_mechanistic_records(b: Builder):
    e=[ev(6)]
    for label in ["INT0","TS1","INT1","INT2","INT2'","INT2''","TS2","TS2'","INT3","TS3","INT4"]:
        b.molecule("m_fig2_"+label.replace("'","p"),label,label,e,identity="unresolved",chirality="unresolved",geometry="unresolved",representation="unresolved",
                   note="Computed Pd complex/transition structure in Figure 2; generic E=CO2Me model and metal coordination are not encoded as a unique SMILES.")
    r=b.base_reaction("figure2_computed_profile",e,"proposed_mechanism",["mechanistic_study"],"Figure 2",
                      "Computed reaction profile for Pd(0)-catalyzed formation of cyclopentenes","DFT-supported proposed catalytic pathway")
    r["participants"]=[b.participant(f"node_{i}","m_fig2_"+label.replace("'","p"),"unknown",e) for i,label in enumerate(["INT0","TS1","INT1","INT2","INT2'","INT2''","TS2","TS2'","INT3","TS3","INT4"],1)]
    values=["INT0 −6.5","TS1 −3.2","INT1 −43.4","INT2 −40.4","INT2' −16.8","INT2'' −24.4","TS2 −34.5","TS2' −24.8","INT3 −61.8","TS3 −45.3","INT4 −85.1"]
    r["notes"]=[sf("Relative free energies ΔG (kcal/mol): "+"; ".join(values),{"metric":"relative free energy","energy_unit":"kcal/mol","values":values},e),
                sf("PCM(DMF)-B3LYP-D3/def2-TZVP//PCM(DMF)-B3LYP-D3/def2-SVP", "PCM(DMF)-B3LYP-D3/def2-TZVP//PCM(DMF)-B3LYP-D3/def2-SVP",e)]
    b.classify_none(r,e,"Figure 2 is a computed mechanistic energy profile, not a standalone observed net synthetic experiment.")
    r["review"]["unresolved"]=["Computed structures and transition states are not encoded as concrete SMILES.","Energy values support a proposed pathway but do not establish an experimentally observed elementary-step sequence."]
    ids=[r["reaction_id"]]
    for label in ["2.2-p-syn","2.2-p-anti","2.2-d-syn"]:
        b.molecule("m_scheme4_"+label.replace(".","_").replace("-","_"),label,label,e,identity="unresolved",chirality="unresolved",geometry="unresolved",representation="unresolved",
                   note="Generic Pd-bound conformer with R/E variables; no concrete unique structure encoded.")
    r=b.base_reaction("scheme4_conformational_model",e,"proposed_mechanism",["mechanistic_study"],"Scheme 4",
                      "Conformational effects on product formation by conditions","author-proposed condition-dependent conformational model")
    r["participants"]=[b.participant("p_syn","m_scheme4_2_2_p_syn","unknown",e),b.participant("p_anti","m_scheme4_2_2_p_anti","unknown",e),b.participant("d_syn","m_scheme4_2_2_d_syn","unknown",e)]
    r["notes"]=[sf("2.2-p-syn (from KH, dppf): access to cyclopropane","proximal-bound syn conformer proposed to access cyclopropane",e),
                sf("2.2-p-anti (from KH, dppf): access to cyclopropane","proximal-bound anti conformer proposed to access cyclopropane",e),
                sf("2.2-d-syn (from NaH, dppf): access to cyclopentene","distal-bound syn conformer proposed to access cyclopentene",e)]
    b.classify_none(r,e,"Scheme 4 is an author-proposed conformational rationale rather than a separate observed experiment.")
    r["review"]["unresolved"]=["The conformers contain generic substituents and Pd coordination and are not encoded as unique SMILES."]
    ids.append(r["reaction_id"])
    return ids


def main():
    dataset=json.loads((RUN/"dataset.json").read_text())
    dataset["documents"][0]["title"]=TITLE
    (RUN/"dataset.json").write_text(json.dumps(dataset,ensure_ascii=False,indent=2)+"\n")
    bundle=st.empty(dataset["documents"])
    b=Builder(bundle)
    common_molecules(b)
    p2=add_literature_and_generic(b)
    p3=add_table1(b)
    p4=add_scope_and_scaleup(b)
    p5=add_mechanism_and_chirality(b)
    p6=add_page6_mechanistic_records(b)
    bundle["molecules"]=list(b.molecules.values())
    bundle["reactions"]=b.reactions
    bundle["source_context"]=[
        {"context_id":f"{DOC}:article_metadata","document_id":DOC,"kind":"article_metadata","evidence_ids":[ev(1,"text")],
         "title":TITLE,"doi":"10.1021/acscatal.1c01847","note":"Title and DOI transcribed from page 1; this is source metadata, not an experiment."},
        {"context_id":f"{DOC}:page2_figure1","document_id":DOC,"kind":"non_reaction_figure","evidence_ids":[ev(2)],
         "description":"Figure 1 shows example drugs/natural products (bilobalide, 6'-fluoroaristeromycin, peramivir, jogyamycin, nitriloprostaglandin I2); no current-paper reaction arrow."},
        {"context_id":f"{DOC}:page6_figure2","document_id":DOC,"kind":"computational_mechanistic_context","evidence_ids":[ev(6)],
         "description":"Figure 2 is a DFT reaction profile; values are retained in the proposed-mechanism record and are not treated as synthetic yields or experimental conditions."},
    ]
    (RUN/"submitted.json").write_text(json.dumps(bundle,ensure_ascii=False,indent=2)+"\n")
    raw={
        "document_id":DOC,"page_range":[1,6],"review_mode":"original rendered page images inspected at original detail",
        "article_title":TITLE,"semantic_accuracy_verified":False,
        "page_inventory":{
            "1":["graphical abstract: generic current-work cycloisomerization and cross-coupling/cyclization arrows","abstract and introduction; no concrete row-wise experiment"],
            "2":["Figure 1 non-reaction molecular examples","Scheme 1 Cazes prior-work arrow","Scheme 1 Ma arrows A and B","Scheme 1 two current-work generic arrows"],
            "3":["Table 1 entries 1-20, each retained as a distinct optimization experiment","ligand structure key and table footnotes","discussion text corroborating selected rows"],
            "4":["Table 2 products 3a-3s, each retained as a distinct scope experiment","footnotes a-e","2 mmol scale-up 2r to 3r","RuCl3/NaIO4 dihydroxylation 3r to 6"],
            "5":["Scheme 2 proposed pathways A-D","Scheme 3 transfer-of-chirality experiments a-c","mechanistic prose and computed values discussion"],
            "6":["Figure 2 DFT reaction profile","Scheme 4 proposed conformational model","no new standalone synthetic experiment"],
        },
        "record_counts":{"total":len(b.reactions),"specific_experiments":sum(r["record_type"]=="specific_experiment" for r in b.reactions),
                         "generic_schemes":sum(r["record_type"]=="generic_scheme" for r in b.reactions),
                         "literature_examples":sum(r["record_type"]=="literature_example" for r in b.reactions),
                         "proposed_mechanisms":sum(r["record_type"]=="proposed_mechanism" for r in b.reactions)},
        "limitations":["No explicitly associated SI was supplied.","Most scope structures are retained as source-identified unresolved molecules rather than unchecked concrete SMILES.","Unmarked Table 2 yield bases and selectivity measurement methods remain unreported."]
    }
    (RUN/"raw_model/pages_0001_0006_extraction.json").write_text(json.dumps(raw,ensure_ascii=False,indent=2)+"\n")
    stereo_diagnostics={
        m["molecule_id"]:st.stereo_inspect(m["smiles"]["normalized_value"])
        for m in bundle["molecules"] if m["representation_type"]=="smiles"
    }
    (RUN/"raw_model/pages_0001_0006_stereo_diagnostics.json").write_text(
        json.dumps(stereo_diagnostics,ensure_ascii=False,indent=2)+"\n"
    )
    mention_hints=st.source_mentions(dataset,["2a","3a","4a","5","3r","6"],context=180)
    (RUN/"raw_model/pages_0001_0006_mentions.json").write_text(
        json.dumps(mention_hints,ensure_ascii=False,indent=2)+"\n"
    )
    (RUN/"audit_pages_0001_0006.json").write_text(
        json.dumps(dataset_audit(dataset),ensure_ascii=False,indent=2)+"\n"
    )
    coverage={"document_id":DOC,"semantic_accuracy_verified":False,"pages":[
        {"page":1,"state":"attempted","regions":[
            {"region":"top and graphical abstract","finding":"article title/DOI plus two generic current-work pathways; reused the page-2 generic records to avoid duplicate experiments","evidence_ids":[ev(1),ev(1,"text")]},
            {"region":"abstract and introduction columns","finding":"context and cited prior work; no row-wise current-paper experiment","evidence_ids":[ev(1),ev(1,"text")]}],
         "record_ids":[f"{DOC}:scheme1_current_cross_coupling_generic",f"{DOC}:scheme1_current_cycloisomerization_generic"],"unresolved":["Graphical abstract uses R/E/Ar variables; no concrete identities or yields are assigned."]},
        {"page":2,"state":"attempted","regions":[
            {"region":"Figure 1","finding":"non-reaction examples of cyclopentane-containing drugs/natural products; retained as source_context","evidence_ids":[ev(2),ev(2,"text")]},
            {"region":"Scheme 1 upper Cazes and Ma rows","finding":"three cited literature examples, not current-paper experiments","evidence_ids":[ev(2),ev(2,"text")]},
            {"region":"Scheme 1 lower This work row","finding":"two generic current-work reaction arrows (cross-coupling/cyclization and cycloisomerization)","evidence_ids":[ev(2),ev(2,"text")]},
            {"region":"lower text columns","finding":"discussion introducing Table 1; no additional distinct row","evidence_ids":[ev(2),ev(2,"text")] }],
         "record_ids":p2,"unresolved":["Generic R/E/Ar variables and literature-example product mixtures prevent complete structures and mapping."]},
        {"page":3,"state":"attempted","regions":[
            {"region":"Table 1 header and reaction drawings","finding":"common charges, three possible products, and rt default transcribed","evidence_ids":[ev(3),ev(3,"text")]},
            {"region":"Table 1 rows 1-20","finding":"20 distinct optimization experiments with explicit product/recovery cells retained row-wise","evidence_ids":[ev(3),ev(3,"text")]},
            {"region":"Table 1 footnotes and ligand key","finding":"NMR basis, isolated yield exception, base-loading exception, reflux/66 °C exceptions, and ligand identities transcribed","evidence_ids":[ev(3),ev(3,"text")]},
            {"region":"discussion columns","finding":"corroborates selected row interpretations and identifies competing product 5","evidence_ids":[ev(3),ev(3,"text")] }],
         "record_ids":p3,"unresolved":["Reaction time and detailed procedure are absent; blank cells remain unreported.","Exact atom mapping and complete product stereochemistry are unresolved."]},
        {"page":4,"state":"attempted","regions":[
            {"region":"Table 2 header/general conditions","finding":"common Pd2dba3/dppm/ArI/NaHMDS/DMF conditions and overnight time transcribed","evidence_ids":[ev(4),ev(4,"text")]},
            {"region":"Table 2 blue scope box 3a-3m","finding":"13 distinct aryl/heteroaryl scope experiments","evidence_ids":[ev(4),ev(4,"text")]},
            {"region":"Table 2 red scope box 3n-3s and footnotes","finding":"six substituted-allene/lactone scope experiments; temperature, NMR, mixture and dr annotations retained","evidence_ids":[ev(4),ev(4,"text")]},
            {"region":"Scale up and postfunctionalization arrows","finding":"2 mmol synthesis of 3r and separate dihydroxylation to 6","evidence_ids":[ev(4),ev(4,"text")]},
            {"region":"discussion columns","finding":"reports 6 in 68% yield and 1.4:1 dr; no additional distinct experiment","evidence_ids":[ev(4),ev(4,"text")] }],
         "record_ids":p4,"unresolved":["Unmarked scope yield basis and dr methods are not stated on page 4.","Most scope structures remain unresolved instead of receiving unchecked SMILES."]},
        {"page":5,"state":"attempted","regions":[
            {"region":"Scheme 2","finding":"author-proposed pathways A-D recorded as one proposed-mechanism record with arrow inventory","evidence_ids":[ev(5),ev(5,"text")]},
            {"region":"Scheme 3a-c","finding":"three distinct transfer-of-chirality experiments with input ee, conditions, product ee and E:Z results","evidence_ids":[ev(5),ev(5,"text")]},
            {"region":"mechanistic discussion and computation text","finding":"author interpretation retained as notes/context; not promoted to separate observed experiments","evidence_ids":[ev(5),ev(5,"text")] }],
         "record_ids":p5,"unresolved":["Product absolute configurations and selectivity measurement methods are not reported on page 5.","Scheme 2 entities are generic and not encoded as concrete SMILES."]},
        {"page":6,"state":"attempted","regions":[
            {"region":"Figure 2","finding":"DFT energy profile recorded as computational proposed-mechanism context, not a synthetic experiment","evidence_ids":[ev(6),ev(6,"text")]},
            {"region":"Scheme 4","finding":"author-proposed conformational model recorded separately","evidence_ids":[ev(6),ev(6,"text")]},
            {"region":"lower discussion columns","finding":"discussion of barriers, regioselectivity and stereodivergence; no additional standalone experiment","evidence_ids":[ev(6),ev(6,"text")] }],
         "record_ids":p6,"unresolved":["Computed metal complexes/conformers are generic or coordination-rich and not encoded as unique SMILES."]}
    ]}
    (RUN/"coverage/pages_0001_0006.json").write_text(json.dumps(coverage,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(raw["record_counts"],ensure_ascii=False))


if __name__ == "__main__":
    main()
