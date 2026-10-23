#!/usr/bin/env python3
"""Source-derived v3.1 extraction for ChemCatChem 2021, 13, 763, PDF pages 1-6."""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_a2856b5fccd4e90c"
RUN = ROOT / "runs/current"


def ev(page, image=False):
    return f"{DOC}:page{page}:{'image' if image else 'text'}"


def evidence(*pages):
    result = []
    for page in pages:
        result.extend([ev(page), ev(page, True)])
    return result


def f(raw=None, normalized=None, pages=(), method="original_extraction",
      status=None, derivation=None, alternatives=()):
    return st.fact(raw, normalized, evidence(*pages), method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def missing():
    return st.fact()


def conflict(raw, alternatives):
    eids = []
    for alternative in alternatives:
        eids.extend(alternative["evidence_ids"])
    return st.fact(raw, None, tuple(dict.fromkeys(eids)), status="conflict", alternatives=alternatives)


bundle = json.loads((RUN / "submitted.json").read_text())
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
molecules = {}


def add_molecule(mid, label, pages, name=None, smiles=None, chirality="unresolved",
                 geometry="unresolved", representation=None, stereo_note=None):
    if mid in molecules:
        return mid
    eids = evidence(*pages)
    m = st.new_molecule(mid, DOC, eids)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, pages)
    m["names"]["original"] = f(name or label, name or label, pages)
    m["names"]["english"] = f(name or label, name or label, pages)
    m["cas_numbers"] = missing()
    m["identity_note"] = f(
        "Source identity retained without external CAS enrichment.",
        "Source identity retained without external CAS enrichment.", pages,
    )
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = f(
            f"structure drawing/name for {label}", smiles, pages,
            method="structure_derivation",
            derivation="Candidate connectivity SMILES transcribed from the supplied source drawing/name; parsing does not certify source identity.",
        )
        m["stereochemistry"]["identity_key"] = f(
            None, diag["stereo_graph_key"], pages, method="structure_derivation",
            derivation="RDKit stereo_graph_key calculated from the candidate SMILES; this is not a sample identity.",
        )
        m["stereochemistry"]["atom_centers"] = f(
            f"graph diagnostic for source drawing/name {label}", diag["atom_centers"], pages,
            method="structure_derivation",
            derivation="RDKit potential/assigned tetrahedral-center diagnostic from the candidate connectivity; '?' remains unassigned.",
        )
        m["stereochemistry"]["double_bonds"] = f(
            f"graph diagnostic for source drawing/name {label}", diag["double_bonds"], pages,
            method="structure_derivation",
            derivation="RDKit double-bond diagnostic from the candidate connectivity; ring/cumulene output is not promoted to a source E/Z assignment.",
        )
        m["stereochemistry"]["isomeric_smiles"] = missing()
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing()
    m["stereochemistry"]["chirality_status"] = f(
        stereo_note or chirality, chirality, pages,
        method="original_extraction" if chirality in ("racemic", "specified", "unspecified") else "inference",
        derivation=None if chirality in ("racemic", "specified", "unspecified") else "Conservative status assigned because the supplied source does not establish a complete stereochemical identity.",
    )
    m["stereochemistry"]["geometry_status"] = f(
        geometry, geometry, pages, method="inference",
        derivation="E/Z is not asserted unless explicitly printed; cyclic alkene geometry and uncertain substituent geometry are not converted into an E/Z assignment.",
    )
    rep = representation or ("partial" if smiles and chirality == "racemic" else "unresolved")
    m["stereochemistry"]["representation_status"] = f(
        rep, rep, pages, method="inference",
        derivation="Connectivity may be represented, while racemic/unknown stereochemical composition is not encoded as a single stereoisomer.",
    )
    if chirality == "racemic":
        m["stereochemistry"]["mixture"] = f(
            "racemic (optically inactive substrate/product series)",
            {"type": "racemate", "composition": "not quantified"}, pages,
        )
        m["stereochemistry"]["assignment"] = f(
            "racemic; individual absolute configurations not assigned",
            "racemic; absolute configurations unassigned", pages,
        )
    elif chirality == "specified":
        m["stereochemistry"]["assignment"] = f(stereo_note, stereo_note, pages)
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


def reagent(mid, label, pages=(2,), name=None):
    return add_molecule(mid, label, pages, name=name)


# Primary structures that are legible enough for conservative connectivity candidates.
add_molecule("mol_1a", "rac-1a", (2, 4, 5, 6), "racemic allenic alcohol 1a",
             "CC(C)=C=CC(O)(C)c1ccccc1", "racemic", "not_applicable")
add_molecule("mol_2a", "rac-2a", (2, 5), "racemic bromodihydrofuran 2a",
             "CC1(C)OC(C)(c2ccccc2)C=C1Br", "racemic", "not_applicable")
add_molecule("mol_3a", "rac-3a", (2, 3, 4, 6),
             "2,2,5-trimethyl-3,5-diphenyl-2,5-dihydrofuran (3a)",
             "CC1(C)OC(C)(c2ccccc2)C=C1c1ccccc1", "racemic", "not_applicable")
add_molecule("mol_9", "9", (4, 5), "non-halogenated cycloisomerization product 9",
             "CC1(C)OC(C)(c2ccccc2)C=C1", "racemic", "not_applicable")
molecules["mol_1a"]["stereochemistry"]["axial"] = f(
    "allene drawing has two identical methyl substituents at one terminus",
    {"axis": "allene", "stereogenic": False, "basis": "identical terminal substituents"},
    (2, 4, 5, 6), method="structure_derivation",
    derivation="The allene axis itself is not stereogenic because one terminal carbon bears two identical methyl substituents; the carbinol carbon remains stereogenic and the source labels 1a racemic.",
)

scope_smiles = {
    "3b": "CC1(C)OC(c2ccccc2)C=C1c1ccccc1",
    "3c": "CC1(C)OC(c2ccc(Cl)cc2)C=C1c1ccccc1",
    "3d": "COc1ccc(C2OC(C)(C)C(c3ccccc3)=C2)cc1",
    "3j": "CC1(C)OC(C)(c2ccccc2)C=C1c1cc(C)cc(C)c1",
    "3k": "CC1(C)OC(C)(c2ccccc2)C=C1c1ccccc1C",
    "3l": "CC1(C)OC(C)(c2ccccc2)C=C1c1ccc2c(c1)OCO2",
    "3m": "CC1(C)OC(C)(c2ccccc2)C=C1c1ccc(OC)cc1",
    "3n": "CC1(C)OC(C)(c2ccccc2)C=C1c1ccc(C(F)(F)F)cc1",
    "3o": "CC1(C)OC(C)(c2ccccc2)C=C1c1ccc(Cl)cc1",
    "3q": "CC1(C)OC(C)(c2ccccc2)C=C1c1ccoc1",
    "3r": "CC1(C)OC(C)(c2ccccc2)C=C1C#Cc1ccccc1",
}
for label in ["3b", "3c", "3d", "3e", "3f", "3g", "3h", "3i"]:
    add_molecule(f"mol_{label}", label, (3,), f"dihydrofuran {label}",
                 scope_smiles.get(label), "unresolved" if label in {"3g", "3h", "3i"} else "racemic",
                 "unresolved" if label in {"3f", "3i"} else "not_applicable")
for label in ["3j", "3k", "3l", "3m", "3n", "3o", "3p", "3q"]:
    add_molecule(f"mol_{label}", label, (4,), f"dihydrofuran {label}",
                 scope_smiles.get(label), "racemic", "unresolved" if label == "3p" else "not_applicable")
add_molecule("mol_3r", "3r", (4,), "alkynylated dihydrofuran 3r", scope_smiles["3r"],
             "racemic", "not_applicable")

# Allenol series 1b-i are retained as source labels because compact drawings are not safely encoded here.
for suffix in "bcdefghi":
    add_molecule(f"mol_1{suffix}", f"1{suffix}", (3,), f"allenic alcohol 1{suffix}",
                 None, "unspecified", "unresolved")

# Generic structures used only by the generic scheme records.
for mid, label, name, pages in [
    ("mol_generic_allenol", "1 / allenol", "generic allenic alcohol", (2, 3)),
    ("mol_generic_halodihydrofuran", "halodihydrofuran", "generic halodihydrofuran", (2,)),
    ("mol_generic_boronic_acid", "Ar-B(OH)2 / RB(OH)2", "generic boronic acid", (2, 4)),
    ("mol_generic_aryldihydrofuran", "3", "generic substituted dihydrofuran", (2, 3, 4)),
]:
    add_molecule(mid, label, pages, name)

# Reagents, catalysts, ligands, solvents and work-up agents (no CAS enrichment).
for args in [
    ("mol_cpo", "CPO", (2, 3, 4, 5, 6), "chloroperoxidase from Caldariomyces fumago"),
    ("mol_gox", "GOx", (2, 3, 4, 6), "glucose oxidase from Aspergillus niger"),
    ("mol_catalase", "catalase", (2, 3, 4, 6), "catalase from Micrococcus lysodeikticus"),
    ("mol_nabr", "NaBr", (2, 3, 4, 5, 6), "sodium bromide"),
    ("mol_nacl", "NaCl", (3,), "sodium chloride"),
    ("mol_nai", "NaI", (3,), "sodium iodide"),
    ("mol_glucose", "D-glucose / Glc", (2, 3, 4, 5, 6), "D-(+)-glucose monohydrate"),
    ("mol_brij35", "Brij 35", (2, 6), "Brij 35 (Brij L23)"),
    ("mol_heptane", "n-heptane", (2, 6), "n-heptane"),
    ("mol_citrate_buffer", "citrate buffer", (2, 6), "aqueous citrate buffer, pH 4.5"),
    ("mol_phb", "PhB(OH)2", (2, 3, 6), "phenylboronic acid"),
    ("mol_xylylb", "3,5-(Me)2PhB(OH)2", (4, 5), "3,5-xylylboronic acid"),
    ("mol_otolylb", "o-tolylboronic acid", (4,), "o-tolylboronic acid"),
    ("mol_benzodioxolylb", "benzodioxolylboronic acid", (4,), "benzodioxolylboronic acid"),
    ("mol_panisylb", "p-methoxyphenylboronic acid", (4,), "p-methoxyphenylboronic acid"),
    ("mol_pcf3b", "p-CF3-phenylboronic acid", (4,), "p-(trifluoromethyl)phenylboronic acid"),
    ("mol_pclb", "p-chlorophenylboronic acid", (4,), "p-chlorophenylboronic acid"),
    ("mol_cispropenylb", "cis-propenylboronic acid", (4,), "cis-propenylboronic acid"),
    ("mol_furylb", "2-furylboronic acid", (4,), "2-furylboronic acid"),
    ("mol_phenylacetylene", "phenylacetylene", (4,), "phenylacetylene"),
    ("mol_na2pdcl4", "Na2PdCl4", (2, 3, 4, 6), "sodium tetrachloropalladate"),
    ("mol_pdoac2", "Pd(OAc)2", (2,), "palladium(II) acetate"),
    ("mol_pph3", "PPh3", (2,), "triphenylphosphine"),
    ("mol_xphos", "4 / XPhos", (2,), "XPhos"),
    ("mol_sphos", "5 / SPhos", (2, 3, 4, 6), "SPhos"),
    ("mol_ssphos", "6 / sSPhos", (2,), "sulfonated SPhos"),
    ("mol_ligand7", "7", (2,), "water-soluble phosphine ligand 7"),
    ("mol_ligand8", "8", (2,), "water-soluble phosphine ligand 8"),
    ("mol_k2co3", "K2CO3", (2,), "potassium carbonate"),
    ("mol_k3po4", "K3PO4", (2, 3, 4, 6), "tribasic potassium phosphate"),
    ("mol_goxpd", "GOx/Pd nanobiohybrid", (4, 5), "glucose oxidase/palladium nanobiohybrid"),
    ("mol_acetonitrile", "acetonitrile", (6,), "acetonitrile"),
    ("mol_npentane", "n-pentane", (6,), "n-pentane"),
    ("mol_brine", "brine", (6,), "brine"),
    ("mol_na2so4", "sodium sulfate", (6,), "sodium sulfate"),
    ("mol_silica", "SiO2", (6,), "silica gel"),
    ("mol_diethyl_ether", "diethyl ether", (6,), "diethyl ether"),
]:
    reagent(*args)
molecules["mol_cispropenylb"]["stereochemistry"]["geometry_status"] = f(
    "cis-propenylboronic acid", "specified", (4,))
molecules["mol_cispropenylb"]["stereochemistry"]["double_bonds"] = f(
    "cis-propenylboronic acid", [{"configuration": "cis", "location": "propenyl C=C"}], (4,))
molecules["mol_cispropenylb"]["stereochemistry"]["representation_status"] = f(
    "source geometry retained textually; no concrete structure encoded", "not_encoded", (4,),
    method="inference", derivation="The source specifies cis reagent geometry, but this unresolved reagent record has no stereochemical SMILES.")


def participant(pid, mid, role, pages, stage="reaction", charges=None, ref=None, product_role=None):
    p = st.new_participant(pid, mid, role, evidence(*pages), stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, pages)
    for key, (raw, norm, method, derivation) in (charges or {}).items():
        p["charges"][key] = f(raw, norm, pages, method=method, derivation=derivation)
    if ref:
        p["equivalence_reference"] = f(ref, ref, pages)
    return p


def q(raw, number, unit, method="original_extraction", derivation=None):
    return raw, {"number": number, "unit": unit, "raw_text": raw}, method, derivation


def step(step_id, pages, pids, operation, temperature=None, time=None, atmosphere=None,
         stage="reaction"):
    s = st.new_step(step_id, evidence(*pages), stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, pages)
    if temperature is not None:
        s["temperature"] = temperature if isinstance(temperature, dict) else f(temperature[0], temperature[1], pages)
    if time is not None:
        s["time"] = time if isinstance(time, dict) else f(time[0], time[1], pages)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, pages)
    return s


def outcome(oid, kind, targets, raw, norm, pages, product_role="product",
            yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence(*pages))
    o["product_role"] = f(product_role, product_role, pages)
    o["value"] = f(raw, norm, pages)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, pages)
    if basis:
        o["measurement_basis"] = f(basis, basis, pages)
    if scope:
        o["scope"] = f(scope, scope, pages)
    return o


TAX = {"C-C": 2, "C-O": 3, "C-Cl": 8, "C-Br": 9, "C-I": 10}


def classification(labels, pages, limitation):
    changes = []
    normalized_labels = []
    for i, code in enumerate(labels, 1):
        a, b = code.split("-")[-2:]
        change = {"change_id": f"b{i}", "atoms": [{"atom_map": None, "element": a},
                                                     {"atom_map": None, "element": b}],
                  "before_order": 0, "after_order": 1, "change_type": "formed"}
        changes.append(f(None, change, pages, method="structure_derivation",
                         derivation="Partial net bond-change assignment from the supplied reaction drawing; atom mapping is not claimed."))
        normalized_labels.append({"category_number": TAX[code], "code": code,
                                  "supporting_change_ids": [f"b{i}"]})
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": missing(),
        "bond_changes": changes,
        "labels": f(None, normalized_labels, pages, method="structure_derivation",
                    derivation="Taxonomy labels follow the partial source-drawing bond-change assignment."),
        "coverage": f(None, "partial", pages, method="inference",
                      derivation="No complete atom mapping; only defensible new connections are classified."),
        "limitations": [limitation, "Hydrogen/proton transfers and omitted byproducts are not completely mapped."],
        "model_votes": [],
    }


def reaction(rid, pages, record_type="specific_experiment", categories=None, location=None,
             reaction_name=None, transformation=None, labels=(), limitation="Partial source-drawing classification."):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence(*pages))
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories or [], categories or [], pages)
    r["location"] = f(location or f"PDF page(s) {','.join(map(str, pages))}",
                      {"pages": list(pages), "region": location}, pages)
    if reaction_name:
        r["reaction_name"] = f(reaction_name, reaction_name, pages)
    if transformation:
        r["transformation_type"] = f(None, transformation, pages, method="inference",
                                     derivation="Transformation type assigned from the source reaction name/drawing.")
    r["classification"] = classification(labels, pages, limitation) if labels else {
        "standard_id": st.STANDARD, "mapped_reaction_smiles": missing(), "bond_changes": [],
        "labels": f(None, [], pages, method="inference", derivation="No bond label asserted for this record."),
        "coverage": f(None, "not_applicable", pages, method="inference",
                      derivation="No concrete observed transformation is classified."),
        "limitations": [limitation], "model_votes": [],
    }
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [],
                   "unresolved": ["No explicitly associated supporting information was supplied; omitted details remain unreported."]}
    r["completeness"] = "partial"
    return r


def add(r):
    bundle["reactions"].append(r)
    return r["reaction_id"]


# Scheme 1a is explicitly prior work; Scheme 1b is the current generic concept.
r = reaction("scheme1a_prior_halocyclization", (2,), "literature_example", [], "Scheme 1a",
             "enzymatic allene halocyclization (2017)", "halocyclization", ("C-O",),
             "Generic X substituent prevents a specific C-X taxonomy assignment.")
r["participants"] = [participant("s1a_allenol", "mol_generic_allenol", "reactant", (2,)),
                     participant("s1a_product", "mol_generic_halodihydrofuran", "product", (2,), product_role="product")]
scheme1a_id = add(r)
r = reaction("scheme1b_generic_onepot", (2,), "generic_scheme", [], "Scheme 1b",
             "traceless activation/cross coupling", "one-pot halocyclization/cross coupling", ("C-O", "C-C"))
r["participants"] = [participant("s1b_allenol", "mol_generic_allenol", "reactant", (2,)),
                     participant("s1b_boronic", "mol_generic_boronic_acid", "reactant", (2,)),
                     participant("s1b_product", "mol_generic_aryldihydrofuran", "product", (2,), product_role="product")]
scheme1b_id = add(r)

# Scheme 2 and 3 generic scope arrows.
r = reaction("scheme2_generic_allenol_scope", (3,), "generic_scheme", [], "Scheme 2 top arrow",
             "structural variation of the allenic coupling partner", "one-pot halocyclization/Suzuki coupling", ("C-O", "C-C"))
r["participants"] = [participant("s2_allenol", "mol_generic_allenol", "reactant", (3,)),
                     participant("s2_phb", "mol_phb", "reactant", (3,)),
                     participant("s2_product", "mol_generic_aryldihydrofuran", "product", (3,), product_role="product")]
scheme2_id = add(r)
r = reaction("scheme3_generic_boronic_scope", (4,), "generic_scheme", [], "Scheme 3 top arrow",
             "structural variation of the boronic acid coupling partner", "one-pot halocyclization/Suzuki coupling", ("C-O", "C-C"))
r["participants"] = [participant("s3_1a", "mol_1a", "reactant", (4,)),
                     participant("s3_rb", "mol_generic_boronic_acid", "reactant", (4,)),
                     participant("s3_product", "mol_generic_aryldihydrofuran", "product", (4,), product_role="product")]
scheme3_id = add(r)

# Table 1 optimization rows.
table = [
    (1, "mol_na2pdcl4", 5, "mol_pph3", "mol_k2co3", 10, 80, 68, 19, 3, False, False),
    (2, "mol_na2pdcl4", 5, "mol_xphos", "mol_k2co3", 10, 80, 25, 60, 3, False, False),
    (3, "mol_na2pdcl4", 5, "mol_sphos", "mol_k2co3", 10, 80, 4, 82, 3, False, False),
    (4, "mol_na2pdcl4", 5, "mol_ssphos", "mol_k2co3", 10, 80, 16, 71, 3, False, False),
    (5, "mol_na2pdcl4", 5, "mol_ligand7", "mol_k2co3", 10, 80, 78, 11, 3, False, False),
    (6, "mol_na2pdcl4", 5, "mol_ligand8", "mol_k2co3", 10, 80, 89, 0, 3, False, False),
    (7, "mol_pdoac2", 5, "mol_sphos", "mol_k2co3", 10, 80, 2, 83, 3, False, False),
    (8, "mol_na2pdcl4", 1, "mol_sphos", "mol_k2co3", 10, 80, 7, 82, 18, False, False),
    (9, "mol_na2pdcl4", 5, "mol_sphos", "mol_k2co3", 50, 80, 32, 47, 3, False, False),
    (10, "mol_na2pdcl4", 5, "mol_sphos", "mol_k2co3", 2, 80, 55, 33, 3, False, False),
    (11, "mol_na2pdcl4", 5, "mol_sphos", "mol_k3po4", 10, 80, 0, 88, 3, False, False),
    (12, "mol_na2pdcl4", 5, "mol_sphos", "mol_k3po4", 10, 60, 44, 42, 3, False, False),
    (13, "mol_na2pdcl4", 5, "mol_sphos", "mol_k3po4", 10, 80, 0, 84, 3, True, False),
    (14, "mol_na2pdcl4", 5, "mol_sphos", "mol_k3po4", 10, 80, 0, 85, 3, False, True),
]
for row, pd, pdmol, ligand, base, baseeq, temp, y2a, y3a, hours, no_cat, argon in table:
    labels = ["C-O"] + (["C-Br"] if y2a > 0 else []) + (["C-C"] if y3a > 0 else [])
    r = reaction(f"table1_entry{row:02d}", (2,), categories=["optimization"],
                 location=f"Table 1, entry {row}", reaction_name="optimization of cross coupling parameters",
                 transformation="sequential enzymatic bromocyclization/Suzuki coupling", labels=labels,
                 limitation="Multiple observed products are summarized; no complete atom mapping is claimed.")
    r["general_scheme_ids"] = [scheme1b_id]
    ps = [
        participant("allene", "mol_1a", "reactant", (2,), charges={"amount": q("100 μmol", 100, "umol")}),
        participant("cpo", "mol_cpo", "catalyst", (2,), charges={"volume": q("8.0 μL", 8.0, "uL")}),
        participant("gox", "mol_gox", "catalyst", (2,), charges={"mass": q("2.6 mg", 2.6, "mg")}),
        participant("nabr", "mol_nabr", "reagent", (2,), charges={"amount": q("300 μmol", 300, "umol")}),
        participant("glucose", "mol_glucose", "reagent", (2,), charges={"amount": q("200 μmol", 200, "umol")}),
        participant("brij", "mol_brij35", "additive", (2,), charges={"amount": q("100 μmol", 100, "umol")}),
        participant("heptane", "mol_heptane", "solvent", (2,), charges={"volume": q("5 mL", 5, "mL")}),
        participant("buffer", "mol_citrate_buffer", "solvent", (2,), charges={"volume": q("5 mL", 5, "mL")}),
        participant("phb", "mol_phb", "reactant", (2,), charges={"amount": q("300 μmol", 300, "umol")}),
        participant("pd", pd, "catalyst", (2,), charges={"mol_percent": q(f"{pdmol} mol%", pdmol, "mol%")}, ref="mol_1a"),
        participant("ligand", ligand, "ligand", (2,), charges={"mol_percent": q("15 mol%", 15, "mol%")}, ref="mol_1a"),
        participant("base", base, "reagent", (2,), charges={"equivalents": q(f"{baseeq} eq.", baseeq, "equiv")}, ref="mol_1a"),
        participant("prod2a", "mol_2a", "product", (2,), product_role="residual bromodihydrofuran intermediate"),
        participant("prod3a", "mol_3a", "product", (2,), product_role="desired arylated product"),
    ]
    if not no_cat:
        ps.insert(8, participant("catalase", "mol_catalase", "catalyst", (2,), charges={"volume": q("1.5 μL", 1.5, "uL")}))
    r["participants"] = ps
    r["reaction_scale"] = f("100 μmol allene", {"reference_molecule_id": "mol_1a", "amount": {"number": 100, "unit": "umol"}}, (2,))
    r["solvent_conditions"]["composition"] = f("5 mL heptane + 5 mL citrate buffer (0.1 M, pH 4.5)",
                                                    ["n-heptane", "aqueous citrate buffer"], (2,))
    r["solvent_conditions"]["mixture_ratio"] = f("1:1", "1:1 v/v", (2,))
    r["solvent_conditions"]["concentration"] = f("10 mM Brij 35", {"number": 10, "unit": "mM"}, (2,))
    r["solvent_conditions"]["concentration_basis"] = f("Brij 35 in the emulsion", "Brij 35 in the emulsion", (2,))
    r["steps"] = [step("enzymatic_halocyclization", (2,), ["allene", "cpo", "gox", "nabr", "glucose", "brij", "heptane", "buffer"],
                       "vigorous stirring in aqueous emulsion", ("room temperature", "r.t."),
                       ("2.5 h", {"number": 2.5, "unit": "h"}))]
    if not no_cat:
        r["steps"].append(step("catalase_hold", (2,), ["catalase"], "catalase addition and hold",
                               None, ("15 min", {"number": 15, "unit": "min"})))
    r["steps"].append(step("cross_coupling", (2,), ["phb", "pd", "ligand", "base"], "palladium-based Suzuki coupling",
                           (f"{temp} °C", {"number": temp, "unit": "°C"}),
                           (f"{hours} h", {"number": hours, "unit": "h"}), "argon" if argon else None))
    r["outcomes"] = [
        outcome("yield_2a", "yield", ["mol_2a"], f"{y2a}%", {"number": y2a, "unit": "%"}, (2,),
                "residual intermediate", "NMR yield", "quantitative 1H NMR against dimethyl sulfone",
                "two-step one-pot experiment; bromodihydrofuran remaining after coupling"),
        outcome("yield_3a", "yield", ["mol_3a"], f"{y3a}%", {"number": y3a, "unit": "%"}, (2,),
                "desired product", "NMR yield", "quantitative 1H NMR against dimethyl sulfone",
                "two-step one-pot experiment"),
    ]
    if no_cat:
        r["notes"].append(f("No catalase addition.", "No catalase addition.", (2,)))
    if argon:
        r["notes"].append(f("Argon atmosphere during cross coupling.", "Argon atmosphere during cross coupling.", (2,)))
    add(r)


# NaCl and NaI substitutions discussed as separate experiments.
for salt_mid, salt_label, yld, code in [("mol_nacl", "NaCl", 6, "C-Cl"), ("mol_nai", "NaI", 48, "C-I")]:
    suffix = salt_label.lower()
    r = reaction(f"halide_control_{suffix}", (3,), categories=["optimization", "control"],
                 location=f"Page 3 prose, {salt_label} substitution", reaction_name="halide substitution control",
                 transformation="one-pot halocyclization/Suzuki coupling", labels=("C-O", "C-C"),
                 limitation=f"The transient C-{salt_label[-1]} bond is not a bond in the isolated final product and is not labelled as a net formed bond.")
    r["general_scheme_ids"] = [scheme1b_id]
    r["participants"] = [participant("allene", "mol_1a", "reactant", (3,)),
                         participant("salt", salt_mid, "reagent", (3,)),
                         participant("phb", "mol_phb", "reactant", (3,)),
                         participant("product", "mol_3a", "product", (3,), product_role="desired product")]
    r["steps"] = [step("onepot_sequence", (3,), ["allene", "salt", "phb"],
                       f"optimized protocol with {salt_label} replacing NaBr")]
    r["outcomes"] = [outcome("overall_yield", "yield", ["mol_3a"], f"{yld}%", {"number": yld, "unit": "%"}, (3,),
                               "desired product", "overall two-step yield", None,
                               f"allenol 1a to arylated dihydrofuran 3a with {salt_label}")]
    r["review"]["unresolved"].append("The prose gives only the changed halide and overall yield; exact quantities, times and measurement basis are not restated.")
    add(r)


# Shared scope participants and steps.
scope_yields = {"3a": 85, "3b": 62, "3c": 56, "3d": 60, "3e": 58, "3f": 50, "3g": 88, "3h": 45, "3i": 53}
for label, yld in scope_yields.items():
    pages = (3, 4, 6) if label == "3a" else (3,)
    r = reaction(f"scheme2_{label}", pages, categories=["substrate_scope"],
                 location=f"Scheme 2, product {label}" + ("; Scheme 3 duplicate; representative procedure" if label == "3a" else ""),
                 reaction_name="structural variation of the allenic coupling partner",
                 transformation="one-pot enzymatic bromocyclization/Suzuki coupling", labels=("C-O", "C-C"))
    r["general_scheme_ids"] = [scheme2_id] + ([scheme3_id] if label == "3a" else [])
    substrate = "mol_1a" if label == "3a" else f"mol_1{label[-1]}"
    ps = [participant("allenol", substrate, "reactant", pages),
          participant("phb", "mol_phb", "reactant", pages),
          participant("cpo", "mol_cpo", "catalyst", pages),
          participant("gox", "mol_gox", "catalyst", pages),
          participant("nabr", "mol_nabr", "reagent", pages),
          participant("glucose", "mol_glucose", "reagent", pages),
          participant("catalase", "mol_catalase", "catalyst", pages),
          participant("pd", "mol_na2pdcl4", "catalyst", pages),
          participant("sphos", "mol_sphos", "ligand", pages),
          participant("base", "mol_k3po4", "reagent", pages),
          participant("product", f"mol_{label}", "product", pages, product_role="isolated scope product")]
    if label == "3a":
        charge_map = {
            "allenol": {"mass": q("18.5 mg", 18.5, "mg"), "amount": q("100 μmol", 100, "umol")},
            "phb": {"mass": q("37 mg", 37, "mg"), "amount": q("300 μmol", 300, "umol")},
            "cpo": {"volume": q("8.0 μL", 8, "uL")}, "gox": {"mass": q("2.6 mg", 2.6, "mg")},
            "nabr": {"mass": q("31 mg", 31, "mg"), "amount": q("300 μmol", 300, "umol")},
            "glucose": {"mass": q("39 mg", 39, "mg"), "amount": q("200 μmol", 200, "umol")},
            "catalase": {"volume": q("1.5 μL", 1.5, "uL")},
            "pd": {"mass": q("1.5 mg", 1.5, "mg"), "amount": q("5 μmol", 5, "umol"),
                   "mol_percent": q("5 mol%", 5, "mol%")},
            "sphos": {"mass": q("6.2 mg", 6.2, "mg"), "amount": q("15 μmol", 15, "umol"),
                      "mol_percent": q("15 mol%", 15, "mol%")},
            "base": {"mass": q("209 mg", 209, "mg"), "amount": q("1 mmol", 1, "mmol"),
                     "equivalents": q("10 equiv", 10, "equiv")},
        }
        for p in ps:
            if p["participant_id"] in charge_map:
                for key, (raw, norm, method, derivation) in charge_map[p["participant_id"]].items():
                    p["charges"][key] = f(raw, norm, (6,), method=method, derivation=derivation)
                if any(k in charge_map[p["participant_id"]] for k in ("mol_percent", "equivalents")):
                    p["equivalence_reference"] = f("mol_1a", "mol_1a", (6,))
        ps.extend([
            participant("brij", "mol_brij35", "additive", (6,), charges={"mass": q("123 mg", 123, "mg"), "amount": q("100 μmol", 100, "umol")}),
            participant("heptane", "mol_heptane", "solvent", (6,), charges={"volume": q("5 mL", 5, "mL")}),
            participant("buffer", "mol_citrate_buffer", "solvent", (6,), charges={"volume": q("5 mL", 5, "mL")}),
            participant("acn", "mol_acetonitrile", "workup_agent", (6,), "workup", {"volume": q("6 mL", 6, "mL")}),
            participant("pentane", "mol_npentane", "workup_agent", (6,), "workup",
                        {"volume": q("4 × 20 mL", 80, "mL", "inference", "Four 20 mL extractions sum to 80 mL; the original multiplicity is retained.")}),
            participant("brine", "mol_brine", "workup_agent", (6,), "workup"),
            participant("na2so4", "mol_na2so4", "workup_agent", (6,), "workup"),
            participant("silica", "mol_silica", "purification_agent", (6,), "purification"),
            participant("ether", "mol_diethyl_ether", "purification_agent", (6,), "purification"),
        ])
    r["participants"] = ps
    if label == "3a":
        r["reaction_scale"] = f("allenic alcohol 1a, 100 μmol", {"reference_molecule_id": "mol_1a", "amount": {"number": 100, "unit": "umol"}}, (6,))
        r["solvent_conditions"]["composition"] = f("5 mL n-heptane + 5 mL aqueous citrate buffer (0.1 M, pH 4.5), Brij 35 emulsion",
                                                        ["n-heptane", "aqueous citrate buffer", "Brij 35"], (6,))
        r["solvent_conditions"]["mixture_ratio"] = f("1:1", "1:1 v/v", (6,))
        time1 = conflict(["4 h (Schemes 2/3)", "2.5 h (representative procedure)"], [
            f("4 h", {"number": 4, "unit": "h"}, (3, 4)),
            f("2.5 h", {"number": 2.5, "unit": "h"}, (6,)),
        ])
        time2 = conflict(["16 h (Schemes 2/3)", "3 h (representative procedure)"], [
            f("16 h", {"number": 16, "unit": "h"}, (3, 4)),
            f("3 h", {"number": 3, "unit": "h"}, (6,)),
        ])
        r["steps"] = [
            step("enzymatic_halocyclization", pages, ["allenol", "cpo", "gox", "nabr", "glucose", "brij", "heptane", "buffer"],
                 "vigorous stirring in an aqueous emulsion", ("room temperature", "r.t."), time1, "air"),
            step("catalase_hold", (6,), ["catalase"], "catalase added to decompose excess hydrogen peroxide",
                 None, ("15 min", {"number": 15, "unit": "min"})),
            step("cross_coupling", pages, ["phb", "base", "sphos", "pd"], "heated with vigorous stirring",
                 ("80 °C", {"number": 80, "unit": "°C"}), time2),
            step("phase_separation_and_extraction", (6,), ["acn", "pentane", "brine", "na2so4"],
                 "cool; add acetonitrile; extract with n-pentane; wash with brine; dry over sodium sulfate; concentrate", stage="workup"),
            step("flash_chromatography", (6,), ["silica", "pentane", "ether"],
                 "flash chromatography on SiO2 with n-pentane/diethyl ether 80:1", stage="purification"),
        ]
        r["procedure"]["original_text"] = f(
            "Representative procedure for the one-pot enzymatic halocyclization and palladium-catalyzed cross coupling (page 6).",
            "Representative procedure on page 6; full text retained in document evidence.", (6,))
        r["procedure"]["addition_order"] = f(
            "Brij/buffer + heptane; NaBr, glucose, 1a, CPO; GOx last; catalase; K3PO4, PhB(OH)2, SPhos, Na2PdCl4 last",
            ["emulsion formation", "halocyclization reagents", "GOx initiation", "catalase", "base/boronic acid/ligand/Pd"], (6,))
        r["outcomes"] = [
            outcome("isolated_yield", "yield", ["mol_3a"], "85%", {"number": 85, "unit": "%"}, pages,
                    "desired product", "isolated yield", "isolated after flash chromatography",
                    "overall yield from allenol 1a through the two-step one-pot sequence"),
            outcome("product_mass", "product_amount", ["mol_3a"], "22.2 mg", {"number": 22.2, "unit": "mg"}, (6,),
                    "desired product", None, "weighed isolated product", "clear oil"),
            outcome("product_amount", "product_amount", ["mol_3a"], "84 μmol", {"number": 84, "unit": "umol"}, (6,),
                    "desired product", None, "reported isolated amount", "clear oil"),
        ]
        r["notes"].append(f(
            "Scheme 2/3 print 4 h then 16 h; the representative procedure prints 2.5 h, 15 min catalase hold, then 3 h. Both are retained without choosing.",
            "source condition conflict retained", pages))
        r["review"]["unresolved"].append("It is unclear whether the generic Scheme 2/3 times or the representative page-6 times describe the exact isolated 3a drawing; both remain as alternatives.")
    else:
        r["steps"] = [
            step("enzymatic_halocyclization", (3,), ["allenol", "cpo", "gox", "nabr", "glucose", "catalase"],
                 "scope protocol, stage 1", ("r.t.", "r.t."), ("4 h", {"number": 4, "unit": "h"}), "air"),
            step("cross_coupling", (3,), ["phb", "base", "sphos", "pd"], "scope protocol, stage 2",
                 ("80 °C", {"number": 80, "unit": "°C"}), ("16 h", {"number": 16, "unit": "h"})),
        ]
        r["outcomes"] = [outcome("isolated_yield", "yield", [f"mol_{label}"], f"{yld}%", {"number": yld, "unit": "%"}, (3,),
                                   "desired product", "isolated yield", "isolated (stated for Schemes 2 and 3 in prose)",
                                   "overall yield over the two reaction steps")]
        if label in {"3e", "3f", "3g", "3h", "3i"}:
            r["review"]["unresolved"].append("The compact scope drawing was not converted to a concrete SMILES; label and source image are retained.")
    add(r)


scope3 = [
    ("3j", "mol_xylylb", 84, None), ("3k", "mol_otolylb", 90, None),
    ("3l", "mol_benzodioxolylb", 84, "10 equiv"), ("3m", "mol_panisylb", 72, None),
    ("3n", "mol_pcf3b", 83, None), ("3o", "mol_pclb", 48, None),
    ("3p", "mol_cispropenylb", 29, "15 equiv; 0.3 M in H2O; syringe pump over 15 h"),
    ("3q", "mol_furylb", 76, "15 equiv; 0.3 M in H2O; syringe pump over 15 h"),
]
for label, boronic, yld, special in scope3:
    r = reaction(f"scheme3_{label}", (4,), categories=["substrate_scope"], location=f"Scheme 3, product {label}",
                 reaction_name="structural variation of the boronic acid coupling partner",
                 transformation="one-pot enzymatic bromocyclization/Suzuki coupling", labels=("C-O", "C-C"))
    r["general_scheme_ids"] = [scheme3_id]
    bcharges = {}
    if label == "3l": bcharges["equivalents"] = q("10 eq. RB(OH)2", 10, "equiv")
    if label in {"3p", "3q"}: bcharges["equivalents"] = q("15 eq. RB(OH)2", 15, "equiv")
    r["participants"] = [participant("allenol", "mol_1a", "reactant", (4,)),
                         participant("boronic", boronic, "reactant", (4,), charges=bcharges,
                                     ref="mol_1a" if bcharges else None),
                         participant("cpo", "mol_cpo", "catalyst", (4,)), participant("gox", "mol_gox", "catalyst", (4,)),
                         participant("nabr", "mol_nabr", "reagent", (4,)), participant("glucose", "mol_glucose", "reagent", (4,)),
                         participant("catalase", "mol_catalase", "catalyst", (4,)),
                         participant("pd", "mol_na2pdcl4", "catalyst", (4,)), participant("sphos", "mol_sphos", "ligand", (4,)),
                         participant("base", "mol_k3po4", "reagent", (4,)),
                         participant("product", f"mol_{label}", "product", (4,), product_role="isolated scope product")]
    r["solvent_conditions"]["composition"] = f("aqueous emulsion", "aqueous emulsion", (4,))
    if label in {"3p", "3q"}:
        r["solvent_conditions"]["concentration"] = f("RB(OH)2 (0.3 M in H2O)", {"number": 0.3, "unit": "M"}, (4,))
        r["solvent_conditions"]["concentration_basis"] = f("boronic acid solution in water", "boronic acid solution in water", (4,))
    r["steps"] = [
        step("enzymatic_halocyclization", (4,), ["allenol", "cpo", "gox", "nabr", "glucose", "catalase"],
             "scope protocol, stage 1", ("r.t.", "r.t."), ("4 h", {"number": 4, "unit": "h"}), "air"),
    ]
    if label in {"3p", "3q"}:
        r["steps"].extend([
            step("continuous_boronic_acid_addition", (4,), ["boronic", "base", "sphos", "pd"],
                 "0.3 M aqueous boronic acid solution continuously added by syringe pump during the coupling stage",
                 ("80 °C", {"number": 80, "unit": "°C"}), ("15 h", {"number": 15, "unit": "h"})),
            step("cross_coupling_stage_total", (4,), ["boronic", "base", "sphos", "pd"],
                 "total cross-coupling stage duration; overlaps the 15 h addition and is not an additional sequential 16 h",
                 ("80 °C", {"number": 80, "unit": "°C"}), ("16 h", {"number": 16, "unit": "h"})),
        ])
    else:
        r["steps"].append(step("cross_coupling", (4,), ["boronic", "base", "sphos", "pd"],
                                 "scope protocol, stage 2", ("80 °C", {"number": 80, "unit": "°C"}),
                                 ("16 h", {"number": 16, "unit": "h"})))
    r["outcomes"] = [outcome("isolated_yield", "yield", [f"mol_{label}"], f"{yld}%", {"number": yld, "unit": "%"}, (4,),
                               "desired product", "isolated yield", "isolated (stated for Schemes 2 and 3 in prose)",
                               "overall yield over the two reaction steps")]
    if special:
        r["notes"].append(f(special, special, (4,)))
    if label == "3p":
        r["review"]["unresolved"].append("The reagent is named cis-propenylboronic acid, but the product alkene configuration is not assigned as E/Z in the source text and is not encoded.")
    add(r)


# Scheme 4 concrete Sonogashira sequence.
r = reaction("scheme4_3r", (4,), categories=["substrate_scope", "other"], location="Scheme 4",
             reaction_name="one-pot halogenation/Sonogashira coupling",
             transformation="one-pot enzymatic bromocyclization/Sonogashira coupling", labels=("C-O", "C-C"))
r["general_scheme_ids"] = [scheme1b_id]
r["participants"] = [participant("allenol", "mol_1a", "reactant", (4,)), participant("phenylacetylene", "mol_phenylacetylene", "reactant", (4,)),
                     participant("cpo", "mol_cpo", "catalyst", (4,)), participant("gox", "mol_gox", "catalyst", (4,)),
                     participant("nabr", "mol_nabr", "reagent", (4,)), participant("glucose", "mol_glucose", "reagent", (4,)),
                     participant("pd", "mol_na2pdcl4", "catalyst", (4,)), participant("sphos", "mol_sphos", "ligand", (4,)),
                     participant("base", "mol_k3po4", "reagent", (4,)), participant("product", "mol_3r", "product", (4,), product_role="desired product")]
r["solvent_conditions"]["composition"] = f("aqueous emulsion", "aqueous emulsion", (4,))
r["steps"] = [step("enzymatic_halocyclization", (4,), ["allenol", "cpo", "gox", "nabr", "glucose"], "stage 1",
                   ("r.t.", "r.t."), ("4 h", {"number": 4, "unit": "h"}), "air"),
              step("sonogashira_coupling", (4,), ["phenylacetylene", "pd", "sphos", "base"], "stage 2",
                   ("80 °C", {"number": 80, "unit": "°C"}), ("16 h", {"number": 16, "unit": "h"}))]
r["outcomes"] = [outcome("yield", "yield", ["mol_3r"], "76%", {"number": 76, "unit": "%"}, (4,),
                           "desired product", None, None, "overall one-pot transformation")]
r["review"]["unresolved"].append("Scheme 4 does not state whether 76% is isolated or NMR yield.")
add(r)


# Scheme 5a: one concrete halocyclization gives two products; Figure 1 is characterization only.
r = reaction("scheme5a_nanobiohybrid_halocyclization", (4, 5), categories=["other"], location="Scheme 5a and page-4 prose",
             reaction_name="GOx/Pd nanobiohybrid-catalyzed halocyclization",
             transformation="enzymatic bromocyclization with competing cycloisomerization", labels=("C-O", "C-Br"),
             limitation="C-Br applies to product 2a; the competing product 9 contributes C-O formation only.")
r["participants"] = [participant("allenol", "mol_1a", "reactant", (4, 5)), participant("hybrid", "mol_goxpd", "catalyst", (4, 5)),
                     participant("cpo", "mol_cpo", "catalyst", (4, 5)), participant("nabr", "mol_nabr", "reagent", (4, 5)),
                     participant("glucose", "mol_glucose", "reagent", (4, 5)),
                     participant("prod2a", "mol_2a", "product", (4, 5), product_role="expected bromodihydrofuran"),
                     participant("prod9", "mol_9", "product", (4, 5), product_role="minor non-halogenated cycloisomerization product")]
r["solvent_conditions"]["composition"] = f("aqueous emulsion", "aqueous emulsion", (4, 5))
r["steps"] = [step("nanobiohybrid_halocyclization", (4, 5), ["allenol", "hybrid", "cpo", "nabr", "glucose"],
                   "GOx/Pd nanobiohybrid substituted for GOx in the original halocyclization protocol",
                   ("r.t.", "r.t."))]
r["outcomes"] = [
    outcome("yield_2a", "yield", ["mol_2a"], "69%", {"number": 69, "unit": "%"}, (5,), "expected product"),
    outcome("yield_9", "yield", ["mol_9"], "7%", {"number": 7, "unit": "%"}, (5,), "minor product"),
    outcome("combined_yield", "yield", ["mol_2a", "mol_9"], "overall yield of 76%", {"number": 76, "unit": "%"}, (4,),
            "combined O-heterocyclic products", "combined yield", None, "sum reported in prose"),
]
r["review"]["unresolved"].append("Time, exact charges and yield measurement basis are not printed for Scheme 5a.")
add(r)

# Scheme 5b: separate Suzuki experiment with the same hybrid.
r = reaction("scheme5b_nanobiohybrid_suzuki", (4, 5), categories=["other"], location="Scheme 5b and page-4/5 prose",
             reaction_name="GOx/Pd nanobiohybrid-catalyzed Suzuki coupling",
             transformation="Suzuki coupling", labels=("C-C",))
r["participants"] = [participant("vinyl_bromide", "mol_2a", "reactant", (4, 5)),
                     participant("boronic", "mol_xylylb", "reactant", (4, 5)),
                     participant("hybrid", "mol_goxpd", "catalyst", (4, 5)),
                     participant("product", "mol_3j", "product", (4, 5), product_role="desired arylated product")]
r["solvent_conditions"]["composition"] = f("aqueous emulsion", "aqueous emulsion", (5,))
r["steps"] = [step("nanobiohybrid_suzuki", (4, 5), ["vinyl_bromide", "boronic", "hybrid"],
                   "GOx/Pd nanobiohybrid-catalyzed cross coupling",
                   ("80 °C", {"number": 80, "unit": "°C"}))]
r["outcomes"] = [outcome("yield", "yield", ["mol_3j"], "81%", {"number": 81, "unit": "%"}, (4, 5), "desired product")]
r["review"]["unresolved"].append("Time, exact charges and yield measurement basis are not printed for Scheme 5b.")
add(r)


bundle["source_context"] = [
    {"context_id": f"{DOC}:figure1_nanobiohybrid_characterization", "document_id": DOC,
     "evidence_ids": evidence(4, 5), "kind": "characterization",
     "summary": "Figure 1 contains SEM, XRD, TEM and nanoparticle-size data; Pd content is 40 wt% and particles are about 2 nm. These are not synthetic reaction records."},
    {"context_id": f"{DOC}:scheme4_unproductive_couplings", "document_id": DOC,
     "evidence_ids": evidence(4), "kind": "under_specified_negative_result",
     "summary": "Initial Stille- and Heck-type attempts were unproductive, but no concrete substrate/row/conditions are given; no fabricated experiment record was created."},
    {"context_id": f"{DOC}:nanobiohybrid_true_onepot_attempts", "document_id": DOC,
     "evidence_ids": evidence(5), "kind": "under_specified_experiment_series",
     "summary": "Early true one-pot nanobiohybrid attempts gave an overall 20-30% range, but the number of runs, exact conditions and product-specific basis are not supplied."},
    {"context_id": f"{DOC}:product3a_characterization", "document_id": DOC,
     "evidence_ids": evidence(6), "kind": "product_characterization",
     "summary": "Page 6 reports Rf, 1H/13C NMR, FT-IR and HRMS for isolated 3a; these observations are linked to the preparation but are not separate synthesis experiments."},
    {"context_id": f"{DOC}:si_availability_unassociated", "document_id": DOC,
     "evidence_ids": evidence(1), "kind": "source_limitation",
     "summary": "The article states that supporting information exists online, but no SI file was supplied or explicitly associated in this workspace."},
]

st.write(RUN / "submitted.json", bundle)

stereo_diagnostics = {}
for molecule in bundle["molecules"]:
    smiles = molecule["smiles"]["normalized_value"]
    if smiles:
        stereo_diagnostics[molecule["molecule_id"]] = st.stereo_inspect(smiles)
st.write(RUN / "raw_model/stereo_diagnostics.json", stereo_diagnostics)

inventory = {
    "document_id": DOC,
    "page_range": [1, 6],
    "visual_inspection": "all six supplied original page images inspected",
    "counts": {"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
               "table1_rows": 14, "scheme2_products": 9, "scheme3_new_products": 8},
    "deduplication": ["3a in Scheme 2, Scheme 3 and the representative procedure is one reaction record with all evidence and conflicting times retained."],
    "non_reaction_context": [c["context_id"] for c in bundle["source_context"]],
}
st.write(RUN / "raw_model/page1_6_visual_inventory.json", inventory)

page_records = {i: [] for i in range(1, 7)}
for rxn in bundle["reactions"]:
    cited_pages = {int(e.split(":page", 1)[1].split(":", 1)[0]) for e in rxn["evidence_ids"]}
    for page in cited_pages & set(page_records):
        page_records[page].append(rxn["reaction_id"])

regions = {
    1: [
        {"region": "title/abstract", "finding": "Article title, one-pot halogenation/Suzuki concept and GOx/Pd hybrid summary; no concrete experiment arrow or table on this page.", "evidence_ids": evidence(1)},
        {"region": "introduction", "finding": "Prior literature and study motivation only; SI availability is mentioned but no SI is associated in the workspace.", "evidence_ids": evidence(1)},
    ],
    2: [
        {"region": "Scheme 1a", "finding": "Prior 2017 generic enzymatic allene halocyclization; classified as literature example.", "evidence_ids": evidence(2)},
        {"region": "Scheme 1b", "finding": "Current-work generic one-pot halocyclization/cross-coupling concept.", "evidence_ids": evidence(2)},
        {"region": "Table 1 and footnotes", "finding": "Fourteen distinct optimization entries with row-specific Pd source/loading, ligand, base, temperature, 2a/3a qNMR yields and footnote exceptions.", "evidence_ids": evidence(2)},
    ],
    3: [
        {"region": "discussion of Table 1", "finding": "NaCl and NaI replacement experiments (6% and 48% overall yields) plus interpretation of entries; no separate record for prose repetitions of table rows.", "evidence_ids": evidence(3)},
        {"region": "Scheme 2", "finding": "Generic scope arrow and nine product arrows 3a-i; 3a is deduplicated with page 4 and page 6 evidence.", "evidence_ids": evidence(3)},
        {"region": "stereochemistry prose", "finding": "Chiral substrates were optically inactive and yielded racemic heterocycles; no ee/er/dr values or absolute configurations are reported.", "evidence_ids": evidence(3)},
    ],
    4: [
        {"region": "Scheme 3", "finding": "Generic boronic-acid scope arrow; 3a duplicate plus eight new products 3j-q, including footnoted 10 equiv and 15 equiv syringe-pump cases.", "evidence_ids": evidence(4)},
        {"region": "Scheme 4", "finding": "Concrete one-pot Sonogashira experiment 1a to 3r, 76%; Stille/Heck failures are under-specified source context, not fabricated records.", "evidence_ids": evidence(4)},
        {"region": "nanobiohybrid prose", "finding": "Preparation/characterization discussion and prose evidence for the two Scheme 5 experiments; physical characterization is not synthesis.", "evidence_ids": evidence(4)},
    ],
    5: [
        {"region": "Figure 1", "finding": "SEM, XRD, TEM and particle-size characterization only; retained in source_context, not a reaction.", "evidence_ids": evidence(5)},
        {"region": "Scheme 5a", "finding": "GOx/Pd hybrid halocyclization of 1a to 2a (69%) plus minor 9 (7%), combined 76% in prose.", "evidence_ids": evidence(5)},
        {"region": "Scheme 5b", "finding": "GOx/Pd hybrid Suzuki coupling of 2a with 3,5-xylylboronic acid to 3j (81%).", "evidence_ids": evidence(5)},
        {"region": "left-column prose", "finding": "Under-specified true one-pot hybrid attempts gave 20-30%; retained as context because run count/conditions/product basis are absent.", "evidence_ids": evidence(5)},
    ],
    6: [
        {"region": "general remarks", "finding": "Analytical and purification instrumentation context; no independent synthetic experiments.", "evidence_ids": evidence(6)},
        {"region": "representative procedure", "finding": "Detailed 3a one-pot preparation, work-up, chromatography, 22.2 mg/84 μmol/85%; merged with the 3a scope record.", "evidence_ids": evidence(6)},
        {"region": "3a characterization", "finding": "Rf, NMR, IR and HRMS observations retained as source_context, not additional reactions.", "evidence_ids": evidence(6)},
    ],
}
unresolved = {
    1: ["No SI file is supplied despite the article's availability statement."],
    2: ["Ligand 7/8 identities are retained by source labels; their structures are not promoted to externally verified identities."],
    3: ["Compact structures 1b-i and products 3e-i are not all safely encodable from the page image; labels/drawings remain authoritative."],
    4: ["Product 3p geometry is not promoted from the reagent name to an E/Z product assignment; several yield measurement bases are not printed."],
    5: ["Scheme 5 times, quantities and yield measurement bases are unreported; the 20-30% early attempt series is not separable into run-level records."],
    6: ["Representative procedure times conflict with the generic scope arrows (2.5/3 h vs 4/16 h); both alternatives are retained."],
}
coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": []}
for page in range(1, 7):
    coverage["pages"].append({"page": page, "state": "attempted", "regions": regions[page],
                              "record_ids": sorted(page_records[page]), "unresolved": unresolved[page]})
st.write(RUN / "coverage/pages_0001_0006.json", coverage)

print(json.dumps(inventory, ensure_ascii=False, indent=2))
