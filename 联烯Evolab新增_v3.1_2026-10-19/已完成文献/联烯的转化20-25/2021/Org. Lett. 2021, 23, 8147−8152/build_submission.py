#!/usr/bin/env python3
"""Build source-derived v3.1 records for doc_1d619ea14ae771a9.

The script only serializes the decisions documented in runs/workflow.md.  It does
not alter the extraction tools, schemas, taxonomy, or validation rules.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_1d619ea14ae771a9"
OUT = ROOT / "runs/current"


def ev(page: int, image: bool = True) -> str:
    return f"{DOC}:page{page}:{'image' if image else 'text'}"


def evidence(*pages: int) -> list[str]:
    result = []
    for page in pages:
        result.extend([ev(page), ev(page, False)])
    return result


def set_fact(raw, normalized, pages, method="original_extraction", derivation=None,
             status=None):
    return st.fact(raw, normalized, evidence(*pages), method=method,
                   derivation=derivation, status=status)


def missing_fact(pages=()):
    return st.fact(evidence=evidence(*pages))


def derived(normalized, pages, derivation, raw=None):
    return set_fact(raw, normalized, pages, method="structure_derivation",
                    derivation=derivation)


def molecule(mid, label, pages, name=None, smiles=None, note=None,
             chirality=None, geometry=None, representation=None):
    m = st.new_molecule(mid, DOC, evidence(*pages))
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = set_fact(label, label, pages)
    m["names"]["original"] = set_fact(name or label, name or label, pages)
    m["names"]["english"] = set_fact(name or label, name or label, pages)
    m["names"]["abbreviations"] = missing_fact(pages)
    m["names"]["aliases"] = missing_fact(pages)
    m["cas_numbers"] = missing_fact(pages)
    m["identity_note"] = set_fact(
        note or "Source-labelled identity; CAS was not reported in the supplied main article.",
        note or "Source-labelled identity; CAS was not reported in the supplied main article.",
        pages,
    )
    m["external_properties"] = []
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            smiles, pages,
            "Candidate connectivity transcribed from the visible source drawing/name; RDKit parsing is a graph diagnostic, not source verification.",
            raw={"source_label": label, "candidate_smiles": smiles},
        )
        m["generic_template"] = None
        m["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], pages,
            "RDKit canonical isomeric SMILES for the candidate graph; it does not infer unreported stereochemistry.")
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], pages,
            "RDKit stereo graph key for comparison only; not a sample identity.")
        potential = [p for p in diag["potential_stereo"] if p["specified"] == "Unspecified"]
        if potential:
            chirality = chirality or "unspecified"
            geometry = geometry or "unspecified"
            representation = representation or "partial"
            m["stereochemistry"]["double_bonds"] = derived(
                diag["double_bonds"], pages,
                "Toolkit reports potential unspecified bond stereochemistry; the article supplies no E/Z descriptor for this candidate graph.")
        else:
            chirality = chirality or "achiral"
            geometry = geometry or "not_applicable"
            representation = representation or "encoded"
            m["stereochemistry"]["assignment"] = derived(
                {"basis": "candidate graph has no assigned or potential stereogenic atom/axis in toolkit diagnostics"},
                pages,
                "Graph-based absence of toolkit-detected stereogenic elements; limited for axial and ring stereochemistry.")
        m["stereochemistry"]["atom_centers"] = derived(
            diag["atom_centers"], pages, "RDKit graph diagnostic on the candidate SMILES.")
        if m["stereochemistry"]["double_bonds"]["normalized_value"] is None:
            m["stereochemistry"]["double_bonds"] = derived(
                diag["double_bonds"], pages, "RDKit graph diagnostic on the candidate SMILES.")
        m["stereochemistry"]["axial"] = missing_fact(pages)
        m["stereochemistry"]["mixture"] = missing_fact(pages)
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing_fact(pages)
        m["generic_template"] = None
        chirality = chirality or "unspecified"
        geometry = geometry or "unspecified"
        representation = representation or "not_encoded"
        m["stereochemistry"]["atom_centers"] = missing_fact(pages)
        m["stereochemistry"]["double_bonds"] = missing_fact(pages)
        m["stereochemistry"]["axial"] = missing_fact(pages)
        m["stereochemistry"]["mixture"] = missing_fact(pages)
        m["stereochemistry"]["assignment"] = set_fact(
            "No R/S, axial, E/Z, or ring-cis/trans descriptor was printed for this source-labelled structure.",
            "No source stereochemical assignment; unresolved structure was not converted to achiral or racemic.", pages)
        m["stereochemistry"]["isomeric_smiles"] = missing_fact(pages)
        m["stereochemistry"]["identity_key"] = missing_fact(pages)
    m["stereochemistry"]["chirality_status"] = derived(
        chirality, pages, "Status assigned conservatively after inspection of the source depiction and available graph diagnostics.")
    m["stereochemistry"]["geometry_status"] = derived(
        geometry, pages, "Status assigned conservatively after inspection of the source depiction and available graph diagnostics.")
    m["stereochemistry"]["representation_status"] = derived(
        representation, pages, "Reports whether source stereochemistry is encoded in the stored candidate representation.")
    return m


def generic_molecule(mid, label, fixed, variable_label, pages, reason):
    m = molecule(mid, label, pages, name=label, note=reason,
                 chirality="unspecified", geometry="unspecified",
                 representation="not_encoded")
    m["representation_type"] = "generic_template"
    m["smiles"] = missing_fact(pages)
    m["generic_template"] = {
        "fixed_structure": set_fact(fixed, {"description": fixed}, pages),
        "variables": [{
            "label": variable_label,
            "attachment_points": set_fact(f"attachment of {variable_label} shown in {label}",
                                           f"attachment of {variable_label} shown in {label}", pages),
            "definition": set_fact(f"source variable {variable_label}", f"source variable {variable_label}", pages),
            "allowed_values": missing_fact(pages),
            "constraints": missing_fact(pages),
        }],
        "instance_mapping": [],
        "reason": set_fact(reason, reason, pages),
    }
    return m


M = {}


def add_m(*args, **kwargs):
    m = molecule(*args, **kwargs)
    M[m["molecule_id"]] = m
    return m["molecule_id"]


def add_generic(*args, **kwargs):
    m = generic_molecule(*args, **kwargs)
    M[m["molecule_id"]] = m
    return m["molecule_id"]


# Core substrates and named reagents/solvents.
add_m("m_allenyl_imide_1", "1", (1, 2, 4), "allenyl imide 1",
      "O=C1OCCN1C(=O)C(C)=C=C",
      note="Candidate connectivity from the visible allenyl-imide drawing. The cumulene graph has unspecified potential bond stereo in toolkit diagnostics; no allene-axis assignment is claimed.",
      chirality="unspecified", geometry="unspecified", representation="partial")
add_m("m_thiobenzamide_2a", "2a", (1, 2, 3, 4), "thiobenzamide 2a", "NC(=S)c1ccccc1")
add_m("m_pbu3", "PBu3", (1, 2, 3, 4), "tributylphosphine", "CCCCP(CCCC)CCCC")
add_m("m_bnp_ph2", "BnPPh2", (1, 2, 4), "benzyldiphenylphosphine", "P(Cc1ccccc1)(c1ccccc1)c1ccccc1")
add_m("m_phpme2", "PhPMe2", (2,), "dimethylphenylphosphine", "P(C)(C)c1ccccc1")
add_m("m_ph2pme", "Ph2PMe", (2,), "methyldiphenylphosphine", "P(C)(c1ccccc1)c1ccccc1")
add_m("m_pph3", "PPh3", (2,), "triphenylphosphine", "P(c1ccccc1)(c1ccccc1)c1ccccc1")
add_m("m_dce", "ClCH2CH2Cl", (2, 3), "1,2-dichloroethane", "ClCCCl")
add_m("m_ch2cl2", "CH2Cl2", (2,), "dichloromethane", "ClCCl")
add_m("m_thf", "THF", (2,), "tetrahydrofuran", "C1CCOC1")
add_m("m_toluene", "toluene", (2, 3), "toluene", "Cc1ccccc1")
add_m("m_benzene", "benzene", (2,), "benzene", "c1ccccc1")
add_m("m_table1_entry11_solvent_conflict", "Table 1 entry 11 solvent (conflict)", (2,),
      "unresolved Table 1 entry 11 solvent",
      note="Table 1 entry 11 visually reads toluene, while the adjacent prose enumerates benzene among the screened solvents. Identity is retained as a source conflict without selecting either solvent.",
      chirality="unresolved", geometry="unresolved", representation="unresolved")
add_m("m_etoh", "EtOH", (2,), "ethanol", "CCO")
add_m("m_cs2co3", "Cs2CO3", (1, 3), "cesium carbonate", "[Cs+].[Cs+].[O-]C([O-])=O")
add_m("m_oxazolidinone", "2-oxazolidinone byproduct", (2, 3), "2-oxazolidinone", "O=C1OCCN1")
add_m("m_pnitrophenol", "4-nitrophenol byproduct", (3,), "4-nitrophenol", "O=[N+]([O-])c1ccc(O)cc1")
add_m("m_cd3od", "CD3OD", (4,), "methanol-d4", "[2H]OC([2H])([2H])[2H]")
add_m("m_me2pph", "Me2PPh", (3,), "dimethylphenylphosphine", "P(C)(C)c1ccccc1")
add_m("m_et3n", "Et3N", (3,), "triethylamine", "CCN(CC)CC")
add_m("m_benzoyl_chloride", "PhCOCl", (3,), "benzoyl chloride", "O=C(Cl)c1ccccc1")
add_m("m_allenoate_6", "6", (3,), "alpha-benzyl-substituted allenoate 6",
      "C=C=C(Cc1ccccc1)C(=O)OCC",
      note="Candidate connectivity from Scheme 4; no allene-axis assignment is printed.",
      chirality="unspecified", geometry="unspecified", representation="partial")

# Scheme 2 thioamides. 2a is already registered; individual source labels for
# 2b-y are descriptive because Scheme 2 prints a generic 2 and labels products.
thioamides = {
    "b": ("p-Me thioarylamide corresponding to 3b", "NC(=S)c1ccc(C)cc1"),
    "c": ("p-OMe thioarylamide corresponding to 3c", "NC(=S)c1ccc(OC)cc1"),
    "d": ("p-OH thioarylamide corresponding to 3d", "NC(=S)c1ccc(O)cc1"),
    "e": ("p-NH2 thioarylamide corresponding to 3e", "NC(=S)c1ccc(N)cc1"),
    "f": ("p-tert-butyl thioarylamide corresponding to 3f", "NC(=S)c1ccc(C(C)(C)C)cc1"),
    "g": ("p-F thioarylamide corresponding to 3g", "NC(=S)c1ccc(F)cc1"),
    "h": ("p-Cl thioarylamide corresponding to 3h", "NC(=S)c1ccc(Cl)cc1"),
    "i": ("p-Br thioarylamide corresponding to 3i", "NC(=S)c1ccc(Br)cc1"),
    "j": ("p-CF3 thioarylamide corresponding to 3j", "NC(=S)c1ccc(C(F)(F)F)cc1"),
    "k": ("p-NO2 thioarylamide corresponding to 3k", "NC(=S)c1ccc([N+](=O)[O-])cc1"),
    "l": ("p-CN thioarylamide corresponding to 3l", "NC(=S)c1ccc(C#N)cc1"),
    "m": ("m-OMe thioarylamide corresponding to 3m", "NC(=S)c1cc(OC)ccc1"),
    "n": ("m-Cl thioarylamide corresponding to 3n", "NC(=S)c1cc(Cl)ccc1"),
    "o": ("m-Br thioarylamide corresponding to 3o", "NC(=S)c1cc(Br)ccc1"),
    "p": ("m-CF3 thioarylamide corresponding to 3p", "NC(=S)c1cc(C(F)(F)F)ccc1"),
    "q": ("o-Cl thioarylamide corresponding to 3q", "NC(=S)c1ccccc1Cl"),
    "r": ("o-Br thioarylamide corresponding to 3r", "NC(=S)c1ccccc1Br"),
    "s": ("meta,para-dimethyl thioarylamide corresponding to 3s", "NC(=S)c1cc(C)c(C)cc1"),
    "t": ("meta,para-dichloro thioarylamide corresponding to 3t", "NC(=S)c1cc(Cl)c(Cl)cc1"),
    "u": ("2-naphthyl thioamide corresponding to 3u", None),
    "v": ("2-furyl thioamide corresponding to 3v", "NC(=S)c1ccco1"),
    "w": ("2-thienyl thioamide corresponding to 3w", "NC(=S)c1cccs1"),
    "x": ("thiourea corresponding to 3x", "NC(N)=S"),
    "y": ("thiopivalamide corresponding to 3y", "NC(=S)C(C)(C)C"),
}
for letter, (name, smi) in thioamides.items():
    add_m(f"m_thioamide_2{letter}", f"2 (substrate for 3{letter})", (2, 3), name, smi,
          note="Individual substrate identity is inferred from the Scheme 2 product substituent and accompanying prose; the figure labels the substrate series generically as 2.")

# Products are retained as source-labelled identities rather than guessed dense
# structures. Stereo is explicitly unknown, never silently achiral/racemic.
for letter in "abcdefghijklmnopqrstuvwxy":
    pages = (1, 2, 4, 6) if letter == "a" else (2, 3)
    name = "4H-1,3-thiazin-4-one " + f"3{letter}"
    add_m(f"m_product_3{letter}", f"3{letter}", pages, name,
          note="Dense product drawing retained by source label. No R/S, axial, E/Z, or ring-cis/trans descriptor/ratio is printed; complete SMILES is intentionally unresolved.")
add_m("m_sideproduct_3aprime", "3a′", (2,), "2H-pyrrol-2-one side product 3a′",
      note="Side product drawing retained by source label; complete connectivity/stereochemistry is unresolved.")

# Scheme 3 alkynoates and thiazolone products.
alkynoates = {
    "a": ("4-nitrophenyl 3-phenylpropiolate 4a", "O=[N+]([O-])c1ccc(OC(=O)C#Cc2ccccc2)cc1"),
    "p": ("4-nitrophenyl 3-(4-methylphenyl)propiolate", "O=[N+]([O-])c1ccc(OC(=O)C#Cc2ccc(C)cc2)cc1"),
    "q": ("4-nitrophenyl 3-(4-fluorophenyl)propiolate", "O=[N+]([O-])c1ccc(OC(=O)C#Cc2ccc(F)cc2)cc1"),
    "r": ("4-nitrophenyl 3-(4-chlorophenyl)propiolate", "O=[N+]([O-])c1ccc(OC(=O)C#Cc2ccc(Cl)cc2)cc1"),
    "s": ("4-nitrophenyl 3-(3-bromophenyl)propiolate", "O=[N+]([O-])c1ccc(OC(=O)C#Cc2cc(Br)ccc2)cc1"),
    "t": ("4-nitrophenyl 3-(2-chlorophenyl)propiolate", "O=[N+]([O-])c1ccc(OC(=O)C#Cc2ccccc2Cl)cc1"),
}
for letter, (name, smi) in alkynoates.items():
    add_m(f"m_alkynoate_4{letter}", f"4{letter}", (1, 3, 4), name, smi)
for letter in "abcdefghijklmnopqrst":
    pages = (1, 3, 4, 6) if letter == "a" else (3,)
    add_m(f"m_product_5{letter}", f"5{letter}", pages, f"5-alkenyl thiazolone 5{letter}",
          note="Scheme 3 depicts one alkenyl arrangement but prints no E/Z descriptor or ratio. Dense product connectivity and alkene geometry are retained as unresolved rather than guessed.")
add_m("m_product_7", "7", (3,), "6,7-dihydro-5H-pyrano[2,3-d]thiazole 7",
      note="Scheme 4 product drawing has no wedge/dash, R/S, dr, or ring-cis/trans assignment; complete structure is unresolved.")
add_m("m_product_8", "8", (3,), "5-arylidene thiazolone 8",
      note="Scheme 4 product drawing has an alkene but no E/Z descriptor or ratio; complete structure is unresolved.")
add_m("m_thiobenzamide_2aD", "2a-D", (4,), "deuterated thiobenzamide 2a-D",
      note="Scheme 5 shows N–D incorporation but no isolated isotopologue composition beyond 10% D incorporation; structure retained by label.")
add_m("m_product_3aD", "3a-D", (4,), "deuterated 3a-D",
      note="Scheme 5 reports 66% deuterium incorporation into a methyl group; exact isotopologue distribution is not expanded into a unique SMILES.")

# Generic structures required for Scheme 1 and the two method families.
add_generic("m_generic_thioamide_2", "generic thioamide 2", "R-C(=S)-NH2", "R", (1, 2, 3, 4),
            "Source uses R/Ar variables across the thioamide scope.")
add_generic("m_generic_alkynoate_4", "generic 4-nitrophenyl 3-arylpropiolate 4",
            "4-NO2-C6H4-O-C(=O)-C#C-Ar", "Ar", (1, 3, 4),
            "Source uses Ar/R2 variables across the alkynoate scope.")
add_generic("m_generic_product_3", "generic 1,3-thiazin-4-one 3",
            "4H-1,3-thiazin-4-one core bearing the thioamide-derived R group", "R", (1, 2, 4),
            "Generic current-work [3+3] product family.")
add_generic("m_generic_product_5", "generic 5-alkenyl thiazolone 5",
            "5-alkenyl thiazolone core bearing thioamide-derived R and alkynoate-derived Ar", "R/Ar", (1, 3, 4),
            "Generic current-work [3+2] product family.")
add_generic("m_prev_electrophile_top", "previous-work methyl ketimine/related 1C,1C-bis-electrophile (top)",
            "generic N-tosyl electrophile drawn in Scheme 1a", "R", (1,),
            "Cited-work structure is generic and is not a current-paper experiment.")
add_generic("m_prev_product_top", "previous-work [4+1] product family (top)",
            "generic five-membered N-tosyl product shown in Scheme 1a", "R", (1,),
            "Cited-work product is generic and retained without expansion.")
add_generic("m_prev_electrophile_bottom", "previous-work enamine/related 1C,1C-bis-electrophile (bottom)",
            "generic N-tosyl electrophile drawn in Scheme 1a", "R/Ar", (1,),
            "Cited-work structure is generic and is not a current-paper experiment.")
add_generic("m_prev_product_bottom", "previous-work [4+1] product family (bottom)",
            "generic five-membered N-tosyl product shown in Scheme 1a", "R/Ar", (1,),
            "Cited-work product is generic and retained without expansion.")


def participant(pid, mid, role, pages, stage="reaction", amount=None,
                volume=None, equiv=None, molpct=None, reference=None,
                product_role=None):
    p = st.new_participant(pid, mid, role, evidence(*pages), stage=stage)
    if amount is not None:
        p["charges"]["amount"] = set_fact(f"{amount} mmol", {"number": amount, "unit": "mmol"}, pages)
    if volume is not None:
        p["charges"]["volume"] = set_fact(f"{volume} mL", {"number": volume, "unit": "mL"}, pages)
    if equiv is not None:
        p["charges"]["equivalents"] = set_fact(f"{equiv} equiv", {"number": equiv, "unit": "equiv"}, pages)
    if molpct is not None:
        p["charges"]["mol_percent"] = set_fact(f"{molpct} mol %", {"number": molpct, "unit": "mol%"}, pages)
    if reference is not None:
        p["equivalence_reference"] = set_fact(reference, reference, pages)
    if product_role is not None:
        p["product_role"] = set_fact(product_role, product_role, pages)
    return p


def step(step_id, pages, pids, operation, temperature=None, time=None):
    s = st.new_step(step_id, evidence(*pages))
    s["participant_ids"] = list(pids)
    s["operation"] = set_fact(operation, operation, pages)
    if temperature is not None:
        raw, norm = temperature
        s["temperature"] = set_fact(raw, norm, pages)
    if time is not None:
        raw, norm = time
        s["time"] = set_fact(raw, norm, pages)
    return s


def outcome(oid, kind, target, pages, raw, norm, basis=None, yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, target, evidence(*pages))
    o["value"] = set_fact(raw, norm, pages)
    if basis:
        o["measurement_basis"] = set_fact(basis, basis, pages)
    if yield_type:
        o["yield_type"] = set_fact(yield_type, yield_type, pages)
    if scope:
        o["scope"] = set_fact(scope, scope, pages)
    if kind in ("yield", "product_amount"):
        o["product_role"] = set_fact("product", "product", pages)
    return o


def classification(pages, formed_pairs=None, coverage="partial", limitation=None,
                   raw_basis=None):
    if coverage == "not_applicable":
        return {
            "standard_id": st.STANDARD,
            "mapped_reaction_smiles": missing_fact(pages),
            "bond_changes": [],
            "labels": derived([], pages, "No defined product transformation was observed; no formation labels assigned."),
            "coverage": set_fact(raw_basis or "no defined transformation", "not_applicable", pages,
                                 method="inference", derivation="Complex-mixture/no-defined-product row; intended chemistry is not classified."),
            "limitations": [limitation or "No defined product identity; intended reaction was not used for classification."],
            "model_votes": [],
        }
    changes = []
    grouped = {}
    for i, pair in enumerate(formed_pairs or [], 1):
        a, b, order = pair
        cid = f"b{i}"
        change = {"change_id": cid, "atoms": [{"atom_map": None, "element": a},
                                                {"atom_map": None, "element": b}],
                  "before_order": 0, "after_order": order, "change_type": "formed"}
        changes.append(derived(change, pages,
                               "Partial net bond change inferred from the visible product scaffold and/or the authors' explicit bond-formation statement; atom mapping is unresolved."))
        lab = st.bond_label(a, b, order)
        grouped.setdefault(lab["code"], dict(lab, supporting_change_ids=[]))["supporting_change_ids"].append(cid)
    labels = list(grouped.values())
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": missing_fact(pages),
        "bond_changes": changes,
        "labels": derived(labels, pages,
                          "Taxonomy labels calculated from the listed partial formed-bond changes; no complete atom mapping is claimed.",
                          raw=raw_basis),
        "coverage": set_fact(raw_basis or "partial source-visible bond changes", coverage, pages,
                             method="inference", derivation="Only source-evident formed bonds are listed; hydrogen balance, leaving groups, and bond-order changes are not atom-completely mapped."),
        "limitations": [limitation or "Partial classification only; complete atom mapping, hydrogen balance, broken bonds, and bond-order changes remain unresolved."],
        "model_votes": [],
    }


def reaction(rid, record_type, pages, categories, location, name, transform,
             participants, steps, outcomes, formed_pairs=None, generic_ids=(),
             notes=(), unresolved=(), class_coverage="partial", class_limit=None,
             raw_class_basis=None, scale=None, solvent=None, solvent_volume=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence(*pages))
    r["record_type"] = record_type
    r["experiment_categories"] = set_fact(categories, list(categories), pages) if categories else missing_fact(pages)
    r["location"] = set_fact(location, {"description": location, "pages": list(pages)}, pages)
    r["reaction_name"] = set_fact(name, name, pages)
    r["transformation_type"] = set_fact(transform, transform, pages)
    r["participants"] = participants
    r["steps"] = steps
    r["outcomes"] = outcomes
    r["general_scheme_ids"] = list(generic_ids)
    r["notes"] = [set_fact(n, n, pages) for n in notes]
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved) or [
                       "Machine-extracted record remains needs_review; unavailable SI prevents recovery of full procedure/workup/purification details."]}
    r["completeness"] = "partial"
    r["classification"] = classification(pages, formed_pairs, class_coverage,
                                           class_limit, raw_class_basis)
    if scale:
        ref, amount = scale
        r["reaction_scale"] = set_fact(f"{ref}: {amount} mmol",
                                        {"reference_molecule_id": ref,
                                         "amount": {"number": amount, "unit": "mmol"}}, pages)
    if solvent:
        r["solvent_conditions"]["composition"] = set_fact(solvent, [solvent], pages)
    if solvent_volume is not None:
        r["solvent_conditions"]["total_volume"] = set_fact(
            f"{solvent_volume} mL", {"number": solvent_volume, "unit": "mL"}, pages)
    return r


R = []
GEN33 = f"{DOC}:generic_3plus3"
GEN32 = f"{DOC}:generic_3plus2"

# Previous-work generic arrows from Scheme 1a.
R.append(reaction(
    "literature_scheme1a_top", "literature_example", (1,), (), "Scheme 1a, previous work, top arrow",
    "[4+1] cycloaddition (previous work, ref 10)", "cited-work generic [4+1] cycloaddition",
    [participant("p1", "m_allenyl_imide_1", "reactant", (1,)),
     participant("p2", "m_prev_electrophile_top", "reactant", (1,)),
     participant("p3", "m_pbu3", "catalyst", (1,)),
     participant("p4", "m_prev_product_top", "product", (1,), product_role="product")],
    [], [], formed_pairs=None,
    notes=("Cited previous work; not counted as a current-paper experiment. Scheme reports yields up to 78% for this generic branch.",),
    unresolved=("Generic R-bearing cited-work structures and complete bond changes are not expanded from the compact Scheme 1 drawing.",),
    class_coverage="partial", class_limit="Generic cited-work structures prevent trustworthy bond mapping.",
    raw_class_basis="Scheme 1a previous-work top branch"))
R[-1]["classification"]["labels"] = derived([], (1,), "No taxonomy label asserted because the generic cited-work structures are unresolved.")
R[-1]["classification"]["coverage"] = set_fact("generic cited-work branch", "unresolved", (1,), method="inference", derivation="Generic structures prevent classification.")
R.append(reaction(
    "literature_scheme1a_bottom", "literature_example", (1,), (), "Scheme 1a, previous work, bottom arrow",
    "[4+1] cycloaddition (previous work, ref 10)", "cited-work generic [4+1] cycloaddition",
    [participant("p1", "m_allenyl_imide_1", "reactant", (1,)),
     participant("p2", "m_prev_electrophile_bottom", "reactant", (1,)),
     participant("p3", "m_pbu3", "catalyst", (1,)),
     participant("p4", "m_prev_product_bottom", "product", (1,), product_role="product")],
    [], [], formed_pairs=None,
    notes=("Cited previous work; not counted as a current-paper experiment. Scheme reports yields up to 51% for this generic branch.",),
    unresolved=("Generic R/Ar-bearing cited-work structures and complete bond changes are not expanded from the compact Scheme 1 drawing.",),
    class_coverage="partial", class_limit="Generic cited-work structures prevent trustworthy bond mapping.",
    raw_class_basis="Scheme 1a previous-work bottom branch"))
R[-1]["classification"]["labels"] = derived([], (1,), "No taxonomy label asserted because the generic cited-work structures are unresolved.")
R[-1]["classification"]["coverage"] = set_fact("generic cited-work branch", "unresolved", (1,), method="inference", derivation="Generic structures prevent classification.")

# Current-work generic branches; Scheme 2/3 details are linked to the same family.
R.append(reaction(
    "generic_3plus3", "generic_scheme", (1, 2), (), "Scheme 1b and Scheme 2 header",
    "[3+3] cycloaddition", "phosphine-catalyzed annulation of allenyl imide and thioamides",
    [participant("p1", "m_allenyl_imide_1", "reactant", (1, 2)),
     participant("p2", "m_generic_thioamide_2", "reactant", (1, 2)),
     participant("p3", "m_bnp_ph2", "catalyst", (1, 2)),
     participant("p4", "m_generic_product_3", "product", (1, 2), product_role="product")],
    [], [], formed_pairs=[("C", "N", 1), ("C", "S", 1)],
    notes=("Generic method family; individual experiments are recorded separately.",),
    unresolved=("Generic R group and complete atom mapping are unresolved.",),
    raw_class_basis="Scheme 5/prose: nitrogen addition and sulfur cyclization form the heterocycle"))
R.append(reaction(
    "generic_3plus2", "generic_scheme", (1, 3), (), "Scheme 1c and Scheme 3 header",
    "[3+2] cycloaddition", "phosphine-catalyzed annulation of 4-nitrophenyl 3-arylpropiolates and thioamides",
    [participant("p1", "m_generic_thioamide_2", "reactant", (1, 3)),
     participant("p2", "m_generic_alkynoate_4", "reactant", (1, 3)),
     participant("p3", "m_pbu3", "catalyst", (1, 3)),
     participant("p4", "m_cs2co3", "reagent", (1, 3)),
     participant("p5", "m_generic_product_5", "product", (1, 3), product_role="product")],
    [], [], formed_pairs=[("C", "N", 1), ("C", "S", 1)],
    notes=("Generic method family; individual experiments are recorded separately.",),
    unresolved=("Generic R/Ar groups, proton source, leaving-group balance, and complete atom mapping are unresolved.",),
    raw_class_basis="Page 3 prose explicitly states C(O)-N and alpha-S-C bond formation"))

# Table 1 optimization entries. Entry 12 also represents Scheme 2 product 3a.
table = [
    (1, "m_pbu3", "m_dce", 1, 49, 24),
    (2, "m_phpme2", "m_dce", 1, 50, 21),
    (3, "m_ph2pme", "m_dce", 1, 71, 23),
    (4, "m_pph3", "m_dce", 2, 78, 16),
    (5, "m_bnp_ph2", "m_dce", 1, 82, "trace"),
    (6, "m_bnp_ph2", "m_dce", 12, 89, "trace"),
    (7, None, "m_dce", 12, "complex mixture", None),
    (8, "m_bnp_ph2", "m_ch2cl2", 12, 91, "trace"),
    (9, "m_bnp_ph2", "m_thf", 12, "complex mixture", None),
    (10, "m_bnp_ph2", "m_toluene", 12, 77, "trace"),
    (11, "m_bnp_ph2", "m_table1_entry11_solvent_conflict", 12, 64, "trace"),
    (12, "m_bnp_ph2", "m_etoh", 12, 99, "table dash (–)"),
]
for entry, catalyst, solvent_mid, hours, y3a, yside in table:
    ps = [participant("p1", "m_allenyl_imide_1", "reactant", (2,), amount=0.10),
          participant("p2", "m_thiobenzamide_2a", "reactant", (2,), amount=0.12),
          participant("p_solvent", solvent_mid, "solvent", (2,), volume=1.0)]
    if catalyst:
        ps.append(participant("p_cat", catalyst, "catalyst", (2,), molpct=20,
                              reference="m_allenyl_imide_1"))
    outs = []
    class_cov = "partial"
    pairs = [("C", "N", 1), ("C", "S", 1)]
    class_limit = "Partial classification covers the defined 3a pathway; competing 3a′ bond changes and complete atom/hydrogen balance are unresolved."
    if y3a == "complex mixture":
        outs.append(outcome(f"o{entry}_mixture", "other", [], (2,), "complex mixture", "complex mixture",
                            scope="Table 1 reaction mixture"))
        class_cov = "not_applicable"
        pairs = None
        class_limit = "No defined product was reported; the intended annulation was not used to assign bond labels."
    else:
        ps += [participant("p3a", "m_product_3a", "product", (2,), product_role="major product"),
               participant("p3aprime", "m_sideproduct_3aprime", "product", (2,), product_role="side product")]
        outs.append(outcome(f"o{entry}_3a", "yield", ["m_product_3a"], (2,),
                            f"{y3a}%", {"number": y3a, "unit": "%"},
                            basis="isolated yield", yield_type="isolated",
                            scope=f"Table 1 entry {entry}, product 3a"))
        if isinstance(yside, int):
            outs.append(outcome(f"o{entry}_3aprime", "yield", ["m_sideproduct_3aprime"], (2,),
                                f"{yside}%", {"number": yside, "unit": "%"},
                                basis="isolated yield", yield_type="isolated",
                                scope=f"Table 1 entry {entry}, side product 3a′"))
        elif yside == "trace":
            outs.append(outcome(f"o{entry}_3aprime", "trace", ["m_sideproduct_3aprime"], (2,),
                                "trace", "trace", basis="isolated yield",
                                scope=f"Table 1 entry {entry}, side product 3a′"))
        else:
            outs.append(outcome(f"o{entry}_3aprime", "other", ["m_sideproduct_3aprime"], (2,),
                                "–", "table dash (–)", basis="isolated yield",
                                scope=f"Table 1 entry {entry}, side product 3a′"))
    pids = [p["participant_id"] for p in ps]
    cats = ("optimization", "substrate_scope") if entry == 12 else ("optimization",)
    pages = (2,)
    rec = reaction(
        f"table1_entry_{entry:02d}", "specific_experiment", pages, cats,
        f"Table 1 entry {entry}" + ("; repeated as Scheme 2 product 3a" if entry == 12 else ""),
        "[3+3] cycloaddition optimization", "allenyl-imide/thioamide annulation",
        ps, [step("reaction", pages, pids, "combine under Table 1 conditions",
                  ("rt", "rt"), (f"{hours} h", {"number": hours, "unit": "h"}))],
        outs, formed_pairs=pairs, generic_ids=(GEN33,),
        notes=(("Scheme 2 3a (99%) is the same EtOH/12 h optimized experiment and is linked here rather than duplicated."
                if entry == 12 else "Table 1 footnote states isolated yields."),),
        unresolved=("SI was not supplied; addition order, atmosphere, workup, and purification are unreported in the provided main article.",
                    "Side-product 3a′ complete bond changes remain unresolved."),
        class_coverage=class_cov, class_limit=class_limit,
        raw_class_basis="[3+3] product 3a; Scheme 5 nitrogen addition and sulfur cyclization" if class_cov != "not_applicable" else "complex mixture",
        scale=("m_allenyl_imide_1", 0.10), solvent=solvent_mid, solvent_volume=1.0)
    if entry == 11:
        rec["solvent_conditions"]["composition"] = st.fact(
            ["Table 1 entry 11: toluene", "adjacent prose: benzene"], None, evidence(2),
            method="original_extraction", status="conflict",
            alternatives=[set_fact("Table 1 entry 11: toluene", ["m_toluene"], (2,)),
                          set_fact("adjacent prose solvent list: benzene", ["m_benzene"], (2,))])
        rec["notes"].append(set_fact(
            "Source conflict: Table 1 entry 11 visually reads toluene, whereas the adjacent prose lists benzene among the screened solvents; neither value was selected.",
            "Source conflict retained without resolution.", (2,)))
        rec["review"]["unresolved"].append(
            "Table 1 entry 11 solvent identity conflicts between the table cell (toluene) and adjacent prose (benzene).")
    R.append(rec)

# Scheme 2 substrate scope, excluding 3a because it is Table 1 entry 12.
scope3 = {
    "b": 92, "c": 81, "d": 72, "e": 95, "f": 97, "g": 84, "h": 88,
    "i": 97, "j": 46, "k": 32, "l": 58, "m": 93, "n": 75, "o": 66,
    "p": 55, "q": 66, "r": 58, "s": 92, "t": 72, "u": 94, "v": 71,
    "w": 67, "x": 67, "y": 83,
}
for letter, yld in scope3.items():
    pages = (2, 3) if letter in "mnopqrstuvwxyz" else (2,)
    mid_sub = f"m_thioamide_2{letter}"
    mid_prod = f"m_product_3{letter}"
    ps = [participant("p1", "m_allenyl_imide_1", "reactant", pages, amount=0.10),
          participant("p2", mid_sub, "reactant", pages, amount=0.12),
          participant("p_cat", "m_bnp_ph2", "catalyst", pages, molpct=20,
                      reference="m_allenyl_imide_1"),
          participant("p_solvent", "m_etoh", "solvent", pages, volume=1.0),
          participant("p_product", mid_prod, "product", pages, product_role="product"),
          participant("p_byproduct", "m_oxazolidinone", "product", pages, product_role="leaving-group byproduct")]
    R.append(reaction(
        f"scheme2_3{letter}", "specific_experiment", pages, ("substrate_scope",),
        f"Scheme 2 product 3{letter}", "[3+3] cycloaddition", "allenyl-imide/thioamide annulation",
        ps, [step("reaction", pages, [p["participant_id"] for p in ps],
                  "Scheme 2 standard conditions", ("rt", "rt"),
                  ("12 h", {"number": 12, "unit": "h"}))],
        [outcome(f"o_3{letter}", "yield", [mid_prod], pages, f"{yld}%",
                 {"number": yld, "unit": "%"}, basis="isolated yield",
                 yield_type="isolated", scope=f"Scheme 2 product 3{letter}")],
        formed_pairs=[("C", "N", 1), ("C", "S", 1)], generic_ids=(GEN33,),
        notes=("Scheme 2 footnote: 1 (0.10 mmol), 2 (0.12 mmol), BnPPh2 (20 mol%), EtOH (1.0 mL); isolated yields.",),
        unresolved=("Individual thioamide label is derived from the displayed product substituent because Scheme 2 labels the substrate series generically as 2.",
                    "SI was not supplied; workup, purification, and complete product SMILES/stereochemistry remain unresolved."),
        class_limit="Partial source-evident C-N/C-S ring closure only; complete atom/hydrogen balance and leaving-group bond changes unresolved.",
        raw_class_basis="Scheme 5 nitrogen addition and sulfur 6-endo-trig cyclization",
        scale=("m_allenyl_imide_1", 0.10), solvent="m_etoh", solvent_volume=1.0))

# Scheme 3 standard scope. Products 5a-o vary thioamide; 5p-t vary alkynoate.
scope5 = {
    "a": 81, "b": 66, "c": 92, "d": 87, "e": 66, "f": 68, "g": 54,
    "h": 46, "i": 73, "j": 49, "k": 71, "l": 92, "m": 61, "n": 79,
    "o": 57, "p": 78, "q": 80, "r": 70, "s": 80, "t": 78,
}
thio_for_5 = {
    "a": "m_thiobenzamide_2a", "b": "m_thioamide_2b", "c": "m_thioamide_2c",
    "d": "m_thioamide_2f", "e": "m_thioamide_2g", "f": "m_thioamide_2h",
    "g": "m_thioamide_2i", "h": "m_thioamide_2j", "i": "m_thioamide_2m",
    "j": "m_thioamide_2p", "k": "m_thioamide_2q", "l": "m_thioamide_2v",
    "m": "m_thioamide_2w", "n": "m_thioamide_2u", "o": "m_thioamide_2y",
}
for letter, yld in scope5.items():
    thio_mid = thio_for_5.get(letter, "m_thiobenzamide_2a")
    alk_mid = f"m_alkynoate_4{letter}" if letter in "pqrst" else "m_alkynoate_4a"
    prod_mid = f"m_product_5{letter}"
    pages = (3,)
    ps = [participant("p2", thio_mid, "reactant", pages, amount=0.12),
          participant("p4", alk_mid, "reactant", pages, amount=0.10),
          participant("p_base", "m_cs2co3", "reagent", pages, amount=0.12, equiv=1.2,
                      reference=alk_mid),
          participant("p_cat", "m_pbu3", "catalyst", pages, molpct=30, reference=alk_mid),
          participant("p_solvent", "m_dce", "solvent", pages, volume=1.0),
          participant("p_product", prod_mid, "product", pages, product_role="product"),
          participant("p_byproduct", "m_pnitrophenol", "product", pages, product_role="byproduct")]
    R.append(reaction(
        f"scheme3_5{letter}", "specific_experiment", pages, ("substrate_scope",),
        f"Scheme 3 product 5{letter}", "[3+2] cycloaddition", "alkynoate/thioamide annulation",
        ps, [step("reaction", pages, [p["participant_id"] for p in ps],
                  "Scheme 3 standard conditions", ("rt", "rt"),
                  ("12 h", {"number": 12, "unit": "h"}))],
        [outcome(f"o_5{letter}", "yield", [prod_mid], pages, f"{yld}%",
                 {"number": yld, "unit": "%"}, basis="isolated yield based on 4",
                 yield_type="isolated", scope=f"Scheme 3 product 5{letter}; based on alkynoate 4")],
        formed_pairs=[("C", "N", 1), ("C", "S", 1)], generic_ids=(GEN32,),
        notes=("Scheme 3 footnote: 2 (0.12 mmol), 4 (0.10 mmol), Cs2CO3 (0.12 mmol), PBu3 (30 mol%), DCE (1.0 mL); isolated yields based on 4.",
               "The article states Cs2CO3 neutralizes the 4-nitrophenol byproduct."),
        unresolved=("SI/Table S1 was not supplied; addition order, workup, purification, exact alkene E/Z assignment, and complete atom mapping remain unresolved.",),
        class_limit="Authors explicitly identify C(O)-N and alpha-S-C formation; C-H/proton source, leaving-group cleavage, and bond-order changes remain atom-unresolved.",
        raw_class_basis="Page 3 prose: C(O)-N and alpha-S-C bond formation",
        scale=(alk_mid, 0.10), solvent="m_dce", solvent_volume=1.0))

# Scheme 3 footnote b is a separate 1.0 mmol scale experiment for 5a.
pages = (3,)
ps = [participant("p2", "m_thiobenzamide_2a", "reactant", pages),
      participant("p4", "m_alkynoate_4a", "reactant", pages, amount=1.0),
      participant("p_base", "m_cs2co3", "reagent", pages, equiv=1.2, reference="m_alkynoate_4a"),
      participant("p_cat", "m_pbu3", "catalyst", pages, molpct=30, reference="m_alkynoate_4a"),
      participant("p_solvent", "m_dce", "solvent", pages),
      participant("p_product", "m_product_5a", "product", pages, product_role="product"),
      participant("p_byproduct", "m_pnitrophenol", "product", pages, product_role="byproduct")]
R.append(reaction(
    "scheme3_5a_scaleup", "specific_experiment", pages, ("scale_up",),
    "Scheme 3 product 5a, footnote b scale-up", "[3+2] cycloaddition", "1.0 mmol alkynoate scale-up",
    ps, [step("reaction", pages, [p["participant_id"] for p in ps],
              "Scheme 3 conditions with 1.0 mmol of 4", ("rt", "rt"),
              ("12 h", {"number": 12, "unit": "h"}))],
    [outcome("o_5a_scale", "yield", ["m_product_5a"], pages, "72% (with 1.0 mmol of 4)",
             {"number": 72, "unit": "%"}, basis="isolated yield based on 4",
             yield_type="isolated", scope="Scheme 3 product 5a; 1.0 mmol of 4a")],
    formed_pairs=[("C", "N", 1), ("C", "S", 1)], generic_ids=(GEN32,),
    notes=("Footnote b explicitly changes the alkynoate scale to 1.0 mmol. Other absolute charges/solvent volume are not reprinted for the scale-up and are not silently multiplied.",),
    unresolved=("Absolute thiobenzamide, Cs2CO3, and solvent charges at 1.0 mmol scale are not explicitly reprinted; only Cs2CO3 1.2 equiv and PBu3 30 mol% are retained from the displayed conditions.",
                "SI workup/purification and product E/Z assignment are unavailable."),
    class_limit="Partial C-N/C-S classification; full atom/hydrogen/leaving-group mapping unresolved.",
    raw_class_basis="Page 3 prose: C(O)-N and alpha-S-C bond formation",
    scale=("m_alkynoate_4a", 1.0), solvent="m_dce"))

# Scheme 4 one-pot derivatizations. Yield basis is not stated in the supplied main text.
ps7 = [participant("p2", "m_thiobenzamide_2a", "reactant", (3,), amount=0.10),
       participant("p4", "m_alkynoate_4a", "reactant", (3,), amount=0.10),
       participant("p_cat1", "m_pbu3", "catalyst", (3,), molpct=20, reference="m_alkynoate_4a"),
       participant("p_base", "m_cs2co3", "reagent", (3,), equiv=1.2, reference="m_alkynoate_4a"),
       participant("p_solvent", "m_toluene", "solvent", (3,), volume=1.0),
       participant("p6", "m_allenoate_6", "reactant", (3,), equiv=1.0, reference="m_alkynoate_4a"),
       participant("p_cat2", "m_me2pph", "catalyst", (3,), molpct=10, reference="m_alkynoate_4a"),
       participant("p7", "m_product_7", "product", (3,), product_role="product")]
R.append(reaction(
    "scheme4a_product7", "specific_experiment", (3,), ("derivatization",),
    "Scheme 4a one-pot sequence to 7", "one-pot [3+2] then [4+2] annulation",
    "sequential phosphine-catalyzed annulations", ps7,
    [step("stage1", (3,), ["p2", "p4", "p_cat1", "p_base", "p_solvent"],
          "generate 5a in situ", ("30 °C", {"number": 30, "unit": "°C"}),
          ("12 h", {"number": 12, "unit": "h"})),
     step("stage2", (3,), ["p6", "p_cat2", "p7"],
          "add allenoate 6 and Me2PPh for [4+2] annulation", None,
          ("5 h", {"number": 5, "unit": "h"}))],
    [outcome("o7", "yield", ["m_product_7"], (3,), "46%", {"number": 46, "unit": "%"},
             scope="one-pot sequence from 2a and 4a to 7")],
    formed_pairs=[("C", "N", 1), ("C", "S", 1), ("C", "C", 1), ("C", "O", 1)],
    generic_ids=(GEN32,),
    notes=("Scheme 4a explicitly shows 2a and 4a at 0.10 mmol each, stage-1 PBu3 20 mol%/Cs2CO3 1.2 equiv in toluene (1.0 mL) at 30 °C for 12 h, then Me2PPh 10 mol% and 6 (1.0 equiv) for 5 h.",),
    unresolved=("Stage-2 temperature is not printed and is left unreported; yield measurement basis, workup, purification, dr/relative configuration, and complete atom mapping are unavailable.",),
    class_limit="Partial family-level classification: initial C-N/C-S annulation and source-visible new C-C/C-O ring connections; complete mapping and bond-order changes unresolved.",
    raw_class_basis="Scheme 4a [3+2] then [4+2] annulation",
    scale=("m_alkynoate_4a", 0.10), solvent="m_toluene", solvent_volume=1.0))

ps8 = [participant("p2", "m_thiobenzamide_2a", "reactant", (3,), amount=0.10),
       participant("p4", "m_alkynoate_4a", "reactant", (3,), amount=0.10),
       participant("p_cat1", "m_pbu3", "catalyst", (3,), molpct=20, reference="m_alkynoate_4a"),
       participant("p_base", "m_cs2co3", "reagent", (3,), equiv=1.2, reference="m_alkynoate_4a"),
       participant("p_solvent", "m_toluene", "solvent", (3,), volume=1.0),
       participant("p_acyl", "m_benzoyl_chloride", "reactant", (3,), equiv=1.2, reference="m_alkynoate_4a"),
       participant("p_et3n", "m_et3n", "reagent", (3,), equiv=1.2, reference="m_alkynoate_4a"),
       participant("p8", "m_product_8", "product", (3,), product_role="product")]
R.append(reaction(
    "scheme4b_product8", "specific_experiment", (3,), ("derivatization",),
    "Scheme 4b one-pot sequence to 8", "one-pot [3+2] annulation then formal beta-acylation",
    "sequential phosphine-catalyzed annulation/acylation", ps8,
    [step("stage1", (3,), ["p2", "p4", "p_cat1", "p_base", "p_solvent"],
          "generate 5a in situ", ("30 °C", {"number": 30, "unit": "°C"}),
          ("12 h", {"number": 12, "unit": "h"})),
     step("stage2", (3,), ["p_acyl", "p_et3n", "p8"],
          "add PhCOCl and Et3N for formal beta-acylation", None,
          ("12 h", {"number": 12, "unit": "h"}))],
    [outcome("o8", "yield", ["m_product_8"], (3,), "72%", {"number": 72, "unit": "%"},
             scope="one-pot sequence from 2a and 4a to 8")],
    formed_pairs=[("C", "N", 1), ("C", "S", 1), ("C", "C", 1)],
    generic_ids=(GEN32,),
    notes=("Scheme 4b explicitly shows 2a and 4a at 0.10 mmol each, stage-1 PBu3 20 mol%/Cs2CO3 1.2 equiv in toluene (1.0 mL) at 30 °C for 12 h, then PhCOCl and Et3N (1.2 equiv each) for 12 h.",),
    unresolved=("Stage-2 temperature is not printed and is left unreported; yield measurement basis, workup, purification, product E/Z assignment, and complete atom mapping are unavailable.",),
    class_limit="Partial family-level classification: initial C-N/C-S annulation plus source-visible beta-acyl C-C formation; complete mapping and bond-order changes unresolved.",
    raw_class_basis="Scheme 4b [3+2] annulation then formal beta-acylation",
    scale=("m_alkynoate_4a", 0.10), solvent="m_toluene", solvent_volume=1.0))

# Scheme 5 deuterium-labelling experiments.
psd1 = [participant("p2a", "m_thiobenzamide_2a", "reactant", (4,)),
        participant("pcat", "m_bnp_ph2", "catalyst", (4,), molpct=20, reference="m_thiobenzamide_2a"),
        participant("psol", "m_cd3od", "solvent", (4,), volume=1.0),
        participant("p2ad", "m_thiobenzamide_2aD", "product", (4,), product_role="isotopically labelled product")]
R.append(reaction(
    "scheme5_deuterium_control", "specific_experiment", (4,), ("control", "mechanistic_study"),
    "Scheme 5 deuterium-labelling control", "H/D exchange of thiobenzamide in CD3OD",
    "deuterium-labelling control", psd1,
    [step("reaction", (4,), [p["participant_id"] for p in psd1],
          "thiobenzamide without allenyl imide in CD3OD", ("rt", "rt"),
          ("2 h", {"number": 2, "unit": "h"}))],
    [outcome("o_dinc_2a", "selectivity", ["m_thiobenzamide_2aD"], (4,),
             "10% D incorporation", {"number": 10, "unit": "%"},
             scope="deuterium incorporation into thiobenzamide 2a-D; isotopologue order/distribution otherwise unreported")],
    formed_pairs=[("N", "H", 1)],
    notes=("The article interprets this low exchange as evidence that product methyl deuteration mainly arises from CD3OD participation in proton transfer rather than preformed 2a-D.",),
    unresolved=("No isolated yield, measurement method, exact N-D isotopologue distribution, or workup is reported in the supplied main article.",),
    class_limit="Isotope exchange is represented as N-H taxonomy formation with concomitant N-H loss; isotope-specific mapping is unresolved.",
    raw_class_basis="Scheme 5 shows N-D in 2a-D and 10% incorporation",
    solvent="m_cd3od", solvent_volume=1.0))

psd2 = [participant("p1", "m_allenyl_imide_1", "reactant", (4,)),
        participant("p2a", "m_thiobenzamide_2a", "reactant", (4,)),
        participant("pcat", "m_bnp_ph2", "catalyst", (4,), molpct=20, reference="m_allenyl_imide_1"),
        participant("psol", "m_cd3od", "solvent", (4,), volume=1.0),
        participant("p3ad", "m_product_3aD", "product", (4,), product_role="deuterated product")]
R.append(reaction(
    "scheme5_deuterium_product", "specific_experiment", (4,), ("mechanistic_study",),
    "Scheme 5 deuterium-labelling annulation", "[3+3] cycloaddition in CD3OD with deuterium incorporation",
    "deuterium-labelling annulation", psd2,
    [step("reaction", (4,), [p["participant_id"] for p in psd2],
          "annulation in CD3OD", ("rt", "rt"),
          ("2 h", {"number": 2, "unit": "h"}))],
    [outcome("o_3ad_yield", "yield", ["m_product_3aD"], (4,), "92%", {"number": 92, "unit": "%"},
             scope="Scheme 5 deuterated product 3a-D"),
     outcome("o_3ad_dinc", "selectivity", ["m_product_3aD"], (4,),
             "66% D incorporation into the methyl group", {"number": 66, "unit": "%"},
             scope="deuterium incorporation at the product methyl group; exact isotopologue distribution unreported")],
    formed_pairs=[("C", "N", 1), ("C", "S", 1), ("C", "H", 1)], generic_ids=(GEN33,),
    notes=("The article attributes the 66% methyl deuteration mainly to CD3OD participation in the 1,5-H transfer step.",),
    unresolved=("Yield measurement basis, D-incorporation measurement method, exact isotopologue distribution, substrate charges, workup, and purification are not reported in the supplied main article.",),
    class_limit="Partial C-N/C-S cyclization plus net C-D (C-H taxonomy) formation; isotope/hydrogen loss and full atom mapping unresolved.",
    raw_class_basis="Scheme 5 product 3a-D and prose on 1,5-H transfer",
    solvent="m_cd3od", solvent_volume=1.0))

# Proposed mechanisms are not experimental records.
R.append(reaction(
    "scheme5_proposed_mechanism", "proposed_mechanism", (3, 4), (),
    "Scheme 5 proposed pathway I-IV to 3", "proposed reaction pathway for 1,3-thiazin-4-one 3",
    "author-proposed phosphonium-ketene/6-endo-trig pathway",
    [participant("p1", "m_allenyl_imide_1", "reactant", (3, 4)),
     participant("p2", "m_generic_thioamide_2", "reactant", (3, 4)),
     participant("pcat", "m_bnp_ph2", "catalyst", (3, 4)),
     participant("p3", "m_generic_product_3", "product", (3, 4), product_role="proposed pathway product")],
    [], [], formed_pairs=[("C", "N", 1), ("C", "S", 1)], generic_ids=(GEN33,),
    notes=("Author-proposed intermediates I-IV/III′ and a disfavored 5-exo-trig path to V are mechanistic drawings, not isolated experiments.",),
    unresolved=("Intermediates are not registered as isolated products; complete elementary-step atom mapping is unresolved.",),
    class_limit="Net family-level bonds only; elementary steps and intermediates are author-proposed.",
    raw_class_basis="Scheme 5 proposed N addition and sulfur 6-endo-trig cyclization"))
R.append(reaction(
    "scheme6_proposed_mechanism", "proposed_mechanism", (4,), (),
    "Scheme 6 proposed pathway VI-XI to 5a", "proposed reaction pathway for 5-alkenyl thiazolones",
    "author-proposed phosphonium-ketene/5-exo-trig pathway",
    [participant("p2", "m_thiobenzamide_2a", "reactant", (4,)),
     participant("p4", "m_alkynoate_4a", "reactant", (4,)),
     participant("pcat", "m_pbu3", "catalyst", (4,)),
     participant("pbase", "m_cs2co3", "reagent", (4,)),
     participant("p5", "m_product_5a", "product", (4,), product_role="proposed pathway product")],
    [], [], formed_pairs=[("C", "N", 1), ("C", "S", 1)], generic_ids=(GEN32,),
    notes=("Author-proposed intermediates VI-XI and a disfavored 6-endo-trig product are mechanistic drawings, not isolated experiments.",
           "Page 4 DFT discussion compares direct, water-assisted, synchronized, and bimolecular 1,2-H transfer pathways; these calculations are source context, not separate synthesis experiments."),
    unresolved=("Complete elementary-step atom mapping and proton-source balance remain unresolved; DFT details reside in unavailable SI.",),
    class_limit="Net C-N/C-S bonds only; elementary steps are proposed and hydrogen/leaving-group balance is unresolved.",
    raw_class_basis="Scheme 6 proposed N addition and sulfur 5-exo-trig cyclization"))


template = json.loads((OUT / "submitted.json").read_text())
template["documents"] = json.loads((OUT / "dataset.json").read_text())["documents"]
template["molecules"] = list(M.values())
template["reactions"] = R
template["source_context"] = [
    {"context_id": f"{DOC}:xray_3a", "document_id": DOC, "kind": "structure_confirmation",
     "subject_molecule_ids": ["m_product_3a"], "evidence_ids": evidence(2, 6),
     "raw_value": "3a assigned by NMR/HRMS and confirmed by X-ray; CCDC 2061479",
     "note": "Identity support only; not an additional synthesis experiment."},
    {"context_id": f"{DOC}:xray_5a", "document_id": DOC, "kind": "structure_confirmation",
     "subject_molecule_ids": ["m_product_5a"], "evidence_ids": evidence(3, 6),
     "raw_value": "ORTEP of 5a; CCDC 2061465",
     "note": "Identity support only; not an additional synthesis experiment."},
    {"context_id": f"{DOC}:dft_proton_transfer", "document_id": DOC, "kind": "computational_mechanism",
     "evidence_ids": evidence(4),
     "raw_value": {"direct_barrier_kcal_mol": 38.8, "water_assisted_stepwise_kcal_mol": [18.0, 8.4],
                   "water_assisted_synchronized_kcal_mol": 28.1},
     "note": "DFT discussion; not a synthetic experiment. Detailed calculations are in SI, which was not supplied."},
    {"context_id": f"{DOC}:bibliography_pages_5_6", "document_id": DOC, "kind": "non_experimental_pages",
     "evidence_ids": evidence(5, 6),
     "raw_value": "Author information, acknowledgments, and references; page 6 also gives CCDC identifiers.",
     "note": "No current-paper reaction arrow/table/procedure was found on pages 5-6."},
]
(OUT / "submitted.json").write_text(json.dumps(template, ensure_ascii=False, indent=2) + "\n")

raw = {
    "document_id": DOC,
    "visual_reviewed_pages": [1, 2, 3, 4, 5, 6],
    "table1": [{"entry": e, "catalyst": c, "solvent": s, "time_h": h,
                "yield_3a": y, "yield_3aprime": yp} for e, c, s, h, y, yp in table],
    "scheme2": [{"product": "3a", "yield_percent": 99, "deduplicated_to": f"{DOC}:table1_entry_12"}] +
               [{"product": f"3{k}", "yield_percent": v} for k, v in scope3.items()],
    "scheme3": [{"product": f"5{k}", "yield_percent": v} for k, v in scope5.items()] +
               [{"product": "5a", "yield_percent": 72, "scale": "1.0 mmol of 4a"}],
    "scheme4": [{"product": "7", "yield_percent": 46}, {"product": "8", "yield_percent": 72}],
    "scheme5_deuterium": [{"product": "2a-D", "D_incorporation_percent": 10},
                           {"product": "3a-D", "yield_percent": 92, "D_incorporation_percent": 66}],
    "notes": ["Raw visual transcription; formal records are in submitted.json.",
              "No SI document was supplied or associated.",
              "Scheme 1 previous-work arrows and Schemes 5-6 proposed mechanisms are not current-paper specific experiments."],
}
(OUT / "raw_model").mkdir(exist_ok=True)
(OUT / "raw_model/extracted_rows.json").write_text(json.dumps(raw, ensure_ascii=False, indent=2) + "\n")
stereo_diagnostics = {}
for mid, mol in M.items():
    smi = mol["smiles"]["normalized_value"]
    if smi:
        stereo_diagnostics[mid] = st.stereo_inspect(smi)
(OUT / "raw_model/stereo_diagnostics.json").write_text(
    json.dumps(stereo_diagnostics, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(M), "reactions": len(R),
                  "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in R)},
                 ensure_ascii=False))
