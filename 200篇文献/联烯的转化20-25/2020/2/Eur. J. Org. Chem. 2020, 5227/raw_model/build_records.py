#!/usr/bin/env python3
"""Build the v3.1 submission for EJOC 2020, 5227, PDF pages 1-6.

This is a transparent source-to-record build script. It intentionally leaves
complex drawn structures unresolved rather than inventing SMILES.
"""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DATASET = json.loads((RUN / "dataset.json").read_text())
DOC = "doc_3bed13a021d538c6"


def ev(page, image=False):
    return f"{DOC}:page{page}:{'image' if image else 'text'}"


def F(raw=None, normalized=None, evidence=(), method="original_extraction",
      status=None, derivation=None):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation)


def q(number, unit, raw=None):
    out = {"number": number, "unit": unit}
    if raw is not None:
        out["raw_text"] = raw
    return out


def rangeq(low, high, unit, raw):
    return {"minimum": low, "maximum": high, "unit": unit, "raw_text": raw}


def conflict(raw, alternatives, evidence):
    f = F(raw, None, evidence)
    f["data_status"] = "conflict"
    f["alternatives"] = alternatives
    return f


def direct(raw, normalized, evidence):
    return F(raw, normalized, evidence)


def derived(raw, normalized, evidence, derivation):
    return F(raw, normalized, evidence, method="structure_derivation",
             derivation=derivation)


def inferred(raw, normalized, evidence, derivation):
    return F(raw, normalized, evidence, method="inference", derivation=derivation)


def make_molecule(mid, label, name, evidence, smiles=None, stereo="unresolved",
                  assignment=None, mixture=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = direct(label, label, evidence)
    m["names"]["original"] = direct(name, name, evidence)
    m["names"]["english"] = direct(name, name, evidence)
    m["cas_numbers"] = F()
    if smiles is not None:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(name, smiles, evidence,
                              "Candidate SMILES derived from the explicitly printed chemical name; graph parsing is not source verification.")
        diag = st.stereo_inspect(smiles)
        m["stereochemistry"]["atom_centers"] = derived(
            name, diag["atom_centers"], evidence, "RDKit standard stereo graph diagnostic on the candidate SMILES.")
        m["stereochemistry"]["double_bonds"] = derived(
            name, diag["double_bonds"], evidence, "RDKit standard stereo graph diagnostic on the candidate SMILES.")
        m["stereochemistry"]["axial"] = derived(
            name, [], evidence, "No axial element asserted for this small named reagent/solvent candidate; toolkit axial support is limited.")
        m["stereochemistry"]["mixture"] = F()
        m["stereochemistry"]["assignment"] = derived(
            name, {"source_assignment": "none", "graph_diagnostic_only": True}, evidence,
            "No potential tetrahedral or geometric stereochemistry was assigned by the graph diagnostic.")
        m["stereochemistry"]["chirality_status"] = derived(
            name, "achiral", evidence, "Candidate graph has no potential tetrahedral stereocenter; this is a structure derivation.")
        m["stereochemistry"]["geometry_status"] = derived(
            name, "not_applicable", evidence, "No independently assignable E/Z element in the candidate graph.")
        m["stereochemistry"]["representation_status"] = derived(
            name, "not_applicable", evidence, "The candidate graph has no stereochemical element requiring encoding.")
        m["stereochemistry"]["isomeric_smiles"] = derived(
            name, diag["canonical_isomeric_smiles"], evidence, "Canonical isomeric SMILES from standard stereo diagnostic.")
        m["stereochemistry"]["identity_key"] = derived(
            name, diag["stereo_graph_key"], evidence, "Stereo graph key from standard stereo diagnostic; not a sample identity.")
        m["identity_note"] = inferred(
            "Named small molecule; candidate graph generated from name", "candidate_structure_needs_review", evidence,
            "The paper provides a name/abbreviation, not a source SMILES; candidate remains machine-derived.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F()
        if stereo in ("partially_specified", "mixture"):
            m["stereochemistry"]["atom_centers"] = inferred(
                "Solid/dashed stereobonds are visible in the source drawing", {"relative_stereobonds_shown": True}, evidence,
                "Visual inspection identifies relative stereobonds, but no absolute R/S labels or verified atom indexing are available.")
            m["stereochemistry"]["assignment"] = inferred(
                assignment or "Relative stereochemistry shown; absolute configuration not assigned", 
                {"relative_configuration": "shown", "absolute_configuration": "unassigned"}, evidence,
                "Source drawings use stereobonds without absolute configurational labels.")
        else:
            m["stereochemistry"]["atom_centers"] = inferred(
                "Complex drawn structure not encoded", [], evidence,
                "Potential tetrahedral centers cannot be indexed safely without a verified structure graph.")
            m["stereochemistry"]["assignment"] = F()
        m["stereochemistry"]["double_bonds"] = inferred(
            "Drawn unsaturation; E/Z applicability not fully resolved", [], evidence,
            "The source does not provide E/Z labels and no verified graph is available.")
        m["stereochemistry"]["axial"] = inferred(
            "Allene axis reviewed where present; no R_a/S_a assignment is reported", 
            {"allene_axis_assignment": "unreported_or_not_applicable_unresolved"}, evidence,
            "Allene-containing substrates are drawn without axial labels; duplicate terminal substituents and exact substituent identities are not encoded here.")
        m["stereochemistry"]["mixture"] = direct(mixture, mixture, evidence) if mixture else F()
        m["stereochemistry"]["chirality_status"] = inferred(
            stereo, stereo, evidence,
            "Status preserves visible relative stereobonds/mixture wording or unresolved stereogenicity without inventing R/S or racemic identity.")
        m["stereochemistry"]["geometry_status"] = inferred(
            "E/Z not assigned in source and applicability unresolved", "unresolved", evidence,
            "No verified structure graph is available to establish E/Z applicability for every double bond.")
        m["stereochemistry"]["representation_status"] = inferred(
            "No verified stereochemical SMILES", "not_encoded", evidence,
            "Complex source drawing retained without an invented concrete SMILES.")
        m["stereochemistry"]["isomeric_smiles"] = F()
        m["stereochemistry"]["identity_key"] = F()
        m["identity_note"] = direct(
            "Identity retained by paper compound label and source drawing; concrete SMILES not asserted.",
            "source_label_and_drawing_only", evidence)
    return m


def make_template(mid, label, description, evidence, variables):
    m = make_molecule(mid, label, description, evidence)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": direct(description, {"source_drawing": description}, evidence),
        "variables": [{
            "label": v,
            "attachment_points": direct(f"{v} attachment shown in source drawing", {"label": v, "position": "as_drawn"}, evidence),
            "definition": direct(f"{v} as printed", f"source-defined variable {v}", evidence),
            "allowed_values": F(),
            "constraints": F(),
        } for v in variables],
        "instance_mapping": [],
        "reason": direct("Generic R-group structure in Scheme 1", "generic_source_structure", evidence),
    }
    m["identity_note"] = direct(description, "generic_template_not_concrete_identity", evidence)
    return m


def participant(pid, mid, role, evidence, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage)
    if product_role:
        p["product_role"] = direct(product_role, product_role, evidence)
    return p


def set_charge(p, field, raw, normalized, evidence, method="original_extraction"):
    p["charges"][field] = F(raw, normalized, evidence, method=method,
                             derivation=("Inherited from the general procedure on PDF page 4." if method == "inherited_general_condition" else None))


def set_classification(r, label_code=None, category_number=None, formed=0,
                       change_type="formed", before=0, after=1, evidence=(), note=""):
    c = r["classification"]
    c["mapped_reaction_smiles"] = F()
    if label_code is None:
        c["bond_changes"] = []
        c["labels"] = inferred("No observed classifiable transformation", [], evidence,
                               note or "No observed product or the source item is only contextual.")
        c["coverage"] = inferred("not applicable or unresolved", "not_applicable", evidence,
                                 note or "No observed transformation to classify.")
        c["limitations"] = [note or "No bond-formation label asserted."]
        return
    changes = []
    for i in range(formed):
        changes.append({
            "change_id": f"b{i+1}",
            "atoms": [{"atom_map": None, "element": label_code.split("-")[-2] if label_code.startswith("DB-") else label_code.split("-")[0]},
                      {"atom_map": None, "element": label_code.split("-")[-1]}],
            "before_order": before,
            "after_order": after,
            "change_type": change_type,
        })
    c["bond_changes"] = [derived(note, ch, evidence,
                                  "Net bond type/count inferred from the source transformation drawing; atom mapping remains unresolved.") for ch in changes]
    c["labels"] = derived(note, [{"category_number": category_number, "code": label_code,
                                   "supporting_change_ids": [x["change_id"] for x in changes]}], evidence,
                           "Taxonomy label follows all asserted partial net bond changes.")
    c["coverage"] = inferred("Partial graph-level assignment", "partial", evidence,
                             "Complex structures were not atom-mapped; only visually unambiguous net bond type/count is retained.")
    c["limitations"] = ["No complete atom mapping or verified concrete reactant/product SMILES.",
                         "Bond-order shifts and hydrogen transfers are unresolved unless separately stated.", note]


def reaction(rid, evidence, record_type="specific_experiment", categories=None,
             location=None, name=None, transformation=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = direct(categories or [], categories or [], evidence)
    r["location"] = direct(location or rid, location or rid, evidence)
    r["reaction_name"] = direct(name, name, evidence) if name else F()
    r["transformation_type"] = inferred(
        transformation, transformation, evidence,
        "Transformation wording normalized from the source title/caption; not an externally verified reaction ontology.") if transformation else F()
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def add_step(r, sid, evidence, participant_ids, operation=None, temperature=None,
             time=None, stage="reaction", method="original_extraction"):
    s = st.new_step(sid, evidence, stage)
    s["participant_ids"] = participant_ids
    if operation is not None:
        s["operation"] = F(operation, operation, evidence, method=method,
                           derivation=("Inherited from page-4 general procedure." if method == "inherited_general_condition" else None))
    if temperature is not None:
        raw, norm = temperature
        s["temperature"] = F(raw, norm, evidence, method=method,
                             derivation=("Inherited from page-4 general procedure." if method == "inherited_general_condition" else None))
    if time is not None:
        raw, norm = time
        s["time"] = F(raw, norm, evidence, method=method,
                      derivation=("Inherited from page-4 general procedure." if method == "inherited_general_condition" else None))
    r["steps"].append(s)
    return s


def add_outcome(r, oid, kind, target, raw, norm, evidence, role=None,
                yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, [target] if target else [], evidence)
    o["value"] = direct(raw, norm, evidence)
    if role:
        o["product_role"] = direct(role, role, evidence)
    if yield_type:
        o["yield_type"] = direct(yield_type, yield_type, evidence)
    if basis:
        o["measurement_basis"] = direct(basis, basis, evidence)
    if scope:
        o["scope"] = direct(scope, scope, evidence)
    r["outcomes"].append(o)
    return o


bundle = st.empty(DATASET["documents"])
molecules = []

# Named small molecules with candidate structures checked by standard stereo.
small = [
    ("dce", "DCE", "1,2-dichloroethane", [ev(2), ev(4)], "ClCCCl"),
    ("phme", "PhMe", "toluene", [ev(2)], "Cc1ccccc1"),
    ("thf", "THF", "tetrahydrofuran", [ev(2)], "C1CCOC1"),
    ("mecn", "CH3CN", "acetonitrile", [ev(2)], "CC#N"),
    ("dcm", "DCM", "dichloromethane", [ev(3)], "ClCCl"),
    ("mcpba", "m-CPBA", "meta-chloroperoxybenzoic acid", [ev(3)], "O=C(OO)c1cccc(Cl)c1"),
    ("fecl3", "FeCl3", "iron(III) chloride", [ev(3)], "[Fe+3].[Cl-].[Cl-].[Cl-]"),
    ("npm", "N-phenylmaleimide", "N-phenylmaleimide", [ev(3)], "O=C1C=CC(=O)N1c1ccccc1"),
    ("et2o", "Et2O", "diethyl ether", [ev(4)], "CCOCC"),
]
for mid, label, name, evidence, smiles in small:
    molecules.append(make_molecule(f"{DOC}:mol:{mid}", label, name, evidence, smiles=smiles))

# Named but structurally unresolved catalysts/agents and purification materials.
for mid, label, name, evidence in [
    ("catA", "A", "gold(I) phosphine complex A (structure drawn in Table 1)", [ev(2, True)]),
    ("catB", "B", "gold(I) phosphine complex B (structure drawn in Table 1)", [ev(2, True)]),
    ("catC", "C", "NHC-gold(I) complex C (structure drawn in Table 1)", [ev(2, True), ev(3, True), ev(4)]),
    ("catD", "D", "NHC-gold(I) complex D (structure drawn in Table 1)", [ev(2, True)]),
    ("ph3pauntf2", "Ph3PAuNTf2", "triphenylphosphinegold(I) bis(trifluoromethanesulfonyl)imide", [ev(2)]),
    ("silica", "silica gel", "silica gel", [ev(4)]),
    ("petether", "petroleum ether", "petroleum ether (composition not specified)", [ev(4)]),
    ("generic_au", "[Au]+", "generic gold(I) catalyst in Schemes 1 and 5", [ev(2, True), ev(4, True)]),
    ("ph3pausbf6", "Ph3PAuSbF6", "triphenylphosphinegold(I) hexafluoroantimonate", [ev(2, True)]),
    ("ptcl2", "PtCl2", "platinum(II) chloride", [ev(2, True)]),
    ("cyjau", "CyJohnphosAuSbF6", "CyJohnphos gold(I) hexafluoroantimonate", [ev(2, True)]),
]:
    molecules.append(make_molecule(f"{DOC}:mol:{mid}", label, name, evidence))

# Substrate and product registries. Complex drawings are not converted to speculative SMILES.
product_names = {
    "a": "2-(3-Tosylpropyl)-4,5,6,11b-tetrahydro-3H-dicyclopenta[a,b]naphthalene",
    "b": "Ethyl 4-(4,5,6,11b-tetrahydro-3H-dicyclopenta[a,b]naphthalen-2-yl)butanoate",
    "c": "3-(4,5,6,11b-Tetrahydro-3H-dicyclopenta[a,b]naphthalen-2-yl)propyl 4-nitrobenzoate",
    "d": "2-(3-Phenoxypropyl)-4,5,6,11b-tetrahydro-3H-dicyclopenta[a,b]naphthalene",
    "e": "4-Methyl-N-phenyl-N-(3-(4,5,6,11b-tetrahydro-3H-dicyclopenta[a,b]naphthalen-2-yl)propyl)benzenesulfonamide",
    "f": "2-(3-(4,5,6,11b-Tetrahydro-3H-dicyclopenta[a,b]naphthalen-2-yl)propyl)isoindoline-1,3-dione",
    "g": "2-Methyl-4,5,6,11b-tetrahydro-3H-dicyclopenta[a,b]naphthalene",
    "h": "2-Decyl-4,5,6,11b-tetrahydro-3H-dicyclopenta[a,b]naphthalene",
    "i": "2-(3-Phenylpropyl)-4,5,6,11b-tetrahydro-3H-dicyclopenta[a,b]naphthalene",
    "j": "2-Cyclohexyl-4,5,6,11b-tetrahydro-3H-dicyclopenta[a,b]naphthalene",
    "k": "4-Bromo-N-(2-decyl-4,5,6,11b-tetrahydro-3H-dicyclopenta[a,b]naphthalen-10-yl)-N-methylbenzenesulfonamide",
    "l": "4-Bromo-N-(2-(3-((4-fluorophenyl)sulfonyl)propyl)-4,5,6,11b-tetrahydro-3H-dicyclopenta[a,b]naphthalen-10-yl)-N-methylbenzenesulfonamide",
    "m": "4-Bromo-N-(2-(3-((4-chlorophenyl)sulfonyl)propyl)-4,5,6,11b-tetrahydro-3H-dicyclopenta[a,b]naphthalen-10-yl)-N-methylbenzenesulfonamide",
    "n": "10-(4-Bromophenyl)-2-(3-tosylpropyl)-4,5,6,11b-tetrahydro-3H-dicyclopenta[a,b]naphthalene",
    "o": "2-Methyl-10-phenyl-4,5,6,11b-tetrahydro-3H-dicyclopenta[a,b]naphthalene",
    "p": "2-Methyl-10-(trifluoromethyl)-4,5,6,11b-tetrahydro-3H-dicyclopenta[a,b]naphthalene",
    "q": "2,8-Dimethyl-4,5,6,11b-tetrahydro-3H-dicyclopenta[a,b]naphthalene",
    "r": "2,9-Dimethyl-4,5,6,11b-tetrahydro-3H-dicyclopenta[a,b]naphthalene (10:1 regioisomer mixture)",
    "s": "2-Methyl-4,5,6,13b-tetrahydro-3H-dicyclopenta[a,b]phenanthrene",
    "t": "7-Methyl-3-tosyl-3,5b,8,9,10,11-hexahydrocyclopenta[7,7a]indeno[5,4-e]indole",
    "u": "2-Methyl-4,5,6,10b-tetrahydro-3H-cyclopenta[7,7a]indeno[4,5-b]thiophene",
    "v": "2,3-Dibutyl-4,5,6,11b-tetrahydro-3H-dicyclopenta[a,b]naphthalene (2.6:1 diastereomer mixture)",
    "w": "2,3,4,4a,5,6,7,12b-Octahydro-1H-benzo[a]cyclopenta[d]fluorene (6.1:1 diastereomer mixture)",
    "x": "1,2-Dimethyl-4,5,6,11b-tetrahydro-3H-dicyclopenta[a,b]naphthalene",
}
char_page = {**{x: 4 for x in "abcdef"}, **{x: 5 for x in "ghijklmnop"}, **{x: 6 for x in "qrstuvwx"}}
for x in "abcdefghijklmnopqrstuvwx":
    sub_e = [ev(3, True), ev(3)] if x != "a" else [ev(2, True), ev(2), ev(3, True), ev(3)]
    molecules.append(make_molecule(f"{DOC}:mol:1{x}", f"1{x}",
                                   f"1,3,9-ene-allene-yne substrate 1{x} (structure as drawn)", sub_e))
    pe = [ev(3, True), ev(3), ev(char_page[x]), ev(char_page[x], True)]
    mix = None
    st_status = "partially_specified"
    if x == "r":
        mix = "inseparable 10:1 mixture of regioisomers"
    if x == "v":
        mix = "inseparable 2.6:1 mixture of diastereoisomers"; st_status = "mixture"
    if x == "w":
        mix = "inseparable 6.1:1 mixture of diastereoisomers"; st_status = "mixture"
    molecules.append(make_molecule(f"{DOC}:mol:2{x}", f"2{x}", product_names[x], pe,
                                   stereo=st_status, mixture=mix))

# Additional products/intermediates in Schemes 3-4.
for mid, label, name, evidence, stereo, mix in [
    ("3e", "3e", "diepoxide 3e", [ev(3), ev(3, True)], "partially_specified", "single diastereomer"),
    ("4e", "4e", "diketone 4e", [ev(3), ev(3, True)], "partially_specified", None),
    ("1y", "1y", "2,6-dimethylphenyl-terminated 1,3,9-ene-allene-yne 1y", [ev(3), ev(3, True)], "unresolved", None),
    ("5", "5", "cyclopentadiene intermediate 5 (not isolable)", [ev(3), ev(3, True)], "unresolved", None),
    ("6", "6", "N-phenylmaleimide [4+2] adduct 6", [ev(3), ev(3, True)], "partially_specified", None),
]:
    molecules.append(make_molecule(f"{DOC}:mol:{mid}", label, name, evidence, stereo=stereo, mixture=mix))

# Generic/literature Scheme 1 molecule templates.
for mid, label, desc, variables in [
    ("s1a_r", "Scheme 1a substrate", "generic 1,3,8-ene-allene-yne from Liu prior work", ["R"]),
    ("s1a_p", "Scheme 1a product", "generic bicyclo[4.3.0]nonadiene prior-work product", ["R"]),
    ("s1b_r", "Scheme 1b substrate", "generic 1,4,9-/1,4,10-allenenyne from authors' prior work", ["R1", "R2", "R3"]),
    ("s1b_p", "Scheme 1b product(s)", "generic bicyclo[3.3.1]nonane prior-work products", ["R1", "R2", "R3", "R4"]),
    ("s1c_r", "Scheme 1c substrate", "generic 1,3-allene-ene from Toste/Iwasawa prior work", ["R1", "R2", "R3", "R4", "R5"]),
    ("s1c_p", "Scheme 1c product", "generic functionalized cyclopentadiene prior-work product", ["R1", "R2", "R3", "R4", "R5"]),
    ("s1d_r", "Scheme 1d substrate", "generic 1,6-ene-arylalkyne from Echavarren prior work", ["R1", "R2", "R3", "R4", "Z"]),
    ("s1d_p", "Scheme 1d product", "generic cyclopenta[b]naphthalene prior-work product", ["R1", "R2", "R3", "R4", "Z"]),
    ("s1e_r", "Scheme 1e substrate", "generic 1,3,9-ene-allene-yne of this work", ["R1", "R2", "R3", "Ar/Het"]),
    ("s1e_p", "Scheme 1e product", "generic fused spirocarbocycle of this work", ["R1", "R2", "R3", "Ar/Het"]),
]:
    molecules.append(make_template(f"{DOC}:mol:{mid}", label, desc, [ev(2, True), ev(2)], variables))

bundle["molecules"] = molecules


def mol(short):
    return f"{DOC}:mol:{short}"


reactions = []

# Scheme 1 arrow inventory: four prior-work transformations and this-work generic proposal.
scheme1 = [
    ("scheme1a_prior", "literature_example", "s1a_r", "s1a_p", "ph3pausbf6",
     "Liu prior work: 6-endo-dig/Nazarov cyclization cascade", "5 mol%"),
    ("scheme1b_prior", "literature_example", "s1b_r", "s1b_p", "generic_au",
     "Authors' prior work: gold(I)-catalyzed cyclopropanation/[5+3] cycloaddition", None),
    ("scheme1c_prior", "literature_example", "s1c_r", "s1c_p", None,
     "Toste and Iwasawa prior work: 1,3-allene-ene cycloisomerization", None),
    ("scheme1d_prior", "literature_example", "s1d_r", "s1d_p", "cyjau",
     "Echavarren prior work: intramolecular [4+2] cycloaddition", "2 mol%"),
    ("scheme1e_generic", "generic_scheme", "s1e_r", "s1e_p", "generic_au",
     "This work: sequential cycloisomerization/formal [4+2] cycloaddition", None),
]
for rid, rtype, rm, pm, cm, title, loading in scheme1:
    e = [ev(2), ev(2, True)]
    r = reaction(rid, e, record_type=rtype, categories=[], location="Scheme 1",
                 name=title, transformation="generic cycloisomerization/cycloaddition")
    r["participants"] = [participant("reactant", mol(rm), "reactant", e),
                         participant("product", mol(pm), "product", e, product_role="illustrated product")]
    if cm:
        cp = participant("catalyst", mol(cm), "catalyst", e)
        if loading:
            set_charge(cp, "mol_percent", loading, q(float(loading.split()[0]), "mol%", loading), e)
            cp["equivalence_reference"] = direct("Scheme template reactant", mol(rm), e)
        r["participants"].append(cp)
    if rid == "scheme1c_prior":
        for pid, mid, load in [("cat_gold", "ph3pausbf6", 2), ("cat_platinum", "ptcl2", 10)]:
            cp = participant(pid, mol(mid), "catalyst", e)
            set_charge(cp, "mol_percent", f"{load} mol%", q(load, "mol%", f"{load} mol%"), e)
            cp["equivalence_reference"] = direct("Scheme template reactant", mol(rm), e)
            r["participants"].append(cp)
        r["notes"].append(direct("Ph3PAuSbF6 (2 mol%) or PtCl2 (10 mol%)", "alternative prior-work catalysts", e))
    r["review"]["unresolved"] = ["Generic/prior-work arrow only; no concrete execution or yield is asserted from this paper."]
    if rid == "scheme1e_generic":
        set_classification(r, "C-C", 2, formed=3, evidence=e,
                           note="Generic drawing indicates sequential carbocyclization with three net C-C ring-forming connections; exact atom mapping unresolved.")
    else:
        r["classification"]["coverage"] = inferred("Prior-work generic drawing", "unresolved", e,
                                                     "Not reclassified as a concrete experiment of this paper.")
    reactions.append(r)

# Table 1 optimization entries. Entry 3 also represents the repeated 2a scope/preparation execution.
opt = [
    (1, "catA", "dce", 86), (2, "catB", "dce", 88), (3, "catC", "dce", 90),
    (4, "catD", "dce", 87), (5, "ph3pauntf2", "dce", None),
    (6, "catC", "phme", 82), (7, "catC", "thf", None), (8, "catC", "mecn", None),
]
for entry, cat, solvent, yld in opt:
    e = [ev(2), ev(2, True)]
    if entry == 3:
        e += [ev(3), ev(3, True), ev(4), ev(4, True)]
    r = reaction(f"table1_entry{entry}", e, categories=["optimization"] + (["substrate_scope"] if entry == 3 else []),
                 location=f"Table 1 entry {entry}" + ("; Scheme 2 product 2a; Experimental Section" if entry == 3 else ""),
                 name="Gold-catalyzed sequential cycloisomerization of 1a",
                 transformation="sequential 1,3-ene-allene cycloisomerization/formal [4+2] cycloaddition")
    r["general_scheme_ids"] = [f"{DOC}:scheme1e_generic"]
    ps = participant("substrate", mol("1a"), "reactant", e)
    if entry == 3:
        ps["charges"]["amount"] = conflict(
            ["Scheme 2: 0.3 mol", "Table 1/general procedure: 0.3 mmol"],
            [direct("0.3 mol", q(0.3, "mol", "0.3 mol"), [ev(3), ev(3, True)]),
             direct("0.3 mmol", q(0.3, "mmol", "0.3 mmol"), [ev(2), ev(4)])], e)
    else:
        set_charge(ps, "amount", "0.3 mmol", q(0.3, "mmol", "0.3 mmol"), e)
    pc = participant("catalyst", mol(cat), "catalyst", e)
    set_charge(pc, "mol_percent", "5 mol%", q(5, "mol%", "5 mol%"), e)
    pc["equivalence_reference"] = direct("relative to 1a", mol("1a"), e)
    pv = participant("solvent", mol(solvent), "solvent", e)
    set_charge(pv, "volume", "3 mL", q(3, "mL", "3 mL"), e)
    r["participants"] = [ps, pc, pv]
    r["reaction_scale"] = (conflict(
        ["Scheme 2: 0.3 mol", "Table 1/general procedure: 0.3 mmol"],
        [direct("0.3 mol", {"reference_molecule_id": mol("1a"), "amount": q(0.3, "mol", "0.3 mol")}, [ev(3), ev(3, True)]),
         direct("0.3 mmol", {"reference_molecule_id": mol("1a"), "amount": q(0.3, "mmol", "0.3 mmol")}, [ev(2), ev(4)])], e)
        if entry == 3 else direct("1a (0.3 mmol)", {"reference_molecule_id": mol("1a"), "amount": q(0.3, "mmol", "0.3 mmol")}, e))
    r["solvent_conditions"]["composition"] = direct(solvent.upper() if solvent != "phme" else "PhMe", [mol(solvent)], e)
    r["solvent_conditions"]["total_volume"] = direct("3 mL", q(3, "mL", "3 mL"), e)
    ids = [p["participant_id"] for p in r["participants"]]
    add_step(r, "reaction", e, ids, "stirred", ("room temperature", "room temperature"),
             ("12 h", q(12, "h", "12 h")))
    if yld is not None:
        pp = participant("product", mol("2a"), "product", e, product_role="desired product")
        r["participants"].append(pp)
        add_outcome(r, "yield_2a", "yield", mol("2a"), f"{yld}%", q(yld, "%", f"{yld}%"),
                    e, role="desired product", yield_type="isolated yield", basis="isolated")
        if entry == 3:
            add_outcome(r, "amount_2a", "product_amount", mol("2a"), "109 mg", q(109, "mg", "109 mg"),
                        [ev(4), ev(4, True)], role="desired product", basis="reported product mass")
        set_classification(r, "C-C", 2, formed=3, evidence=e,
                           note="Product drawing shows the sequential formation of three carbocyclic C-C connections; exact atom identities/bond-order shifts remain unresolved.")
    else:
        basis = "TLC and 1H NMR measurements of the crude reaction mixture"
        add_outcome(r, "no_reaction", "no_reaction", None, "–; no reaction", True, e, basis=basis)
        set_classification(r, evidence=e, note="No reaction was observed; no intended-product bond formation is assigned.")
    r["review"]["unresolved"] = ["Complex substrate/product structures are retained by labels and drawings without verified SMILES or atom mapping."]
    if entry == 3:
        r["review"]["unresolved"].append("Scheme 2 prints 0.3 mol, whereas Table 1 and the general procedure print 0.3 mmol; both are preserved as a conflict.")
    reactions.append(r)

# Scheme 2 substrate scope, excluding 2a because it is the same execution as Table 1 entry 3.
yields = dict(zip("bcdefghijklmnopqrstuvwx", [84,95,88,92,96,85,90,81,82,82,70,71,80,83,62,78,77,74,76,66,63,60,64]))
masses = dict(zip("abcdefghijklmnopqrstuvwx", [109,81,118,90,137,114,57,94,79,71,147,138,143,134,74,54,55,55,60,95,45,61,47,45]))
for x, yld in yields.items():
    cp = char_page[x]
    e = list(dict.fromkeys([ev(3), ev(3, True), ev(cp), ev(cp, True), ev(4)]))
    r = reaction(f"scope_2{x}", e, categories=["substrate_scope"],
                 location=f"Scheme 2, 1{x} to 2{x}; product paragraph on PDF page {cp}",
                 name="Gold-catalyzed sequential cycloisomerization of 1,3,9-ene-allene-yne",
                 transformation="sequential 1,3-ene-allene cycloisomerization/formal [4+2] cycloaddition")
    r["general_scheme_ids"] = [f"{DOC}:scheme1e_generic"]
    ps = participant("substrate", mol(f"1{x}"), "reactant", e)
    ps["charges"]["amount"] = conflict(
        ["Scheme 2: 0.3 mol", "general procedure: 0.3 mmol"],
        [direct("0.3 mol", q(0.3, "mol", "0.3 mol"), [ev(3), ev(3, True)]),
         direct("0.3 mmol", q(0.3, "mmol", "0.3 mmol"), [ev(4)])], e)
    pcata = participant("catalyst", mol("catC"), "catalyst", e)
    set_charge(pcata, "mass", "14 mg", q(14, "mg", "14 mg"), [ev(4)], "inherited_general_condition")
    set_charge(pcata, "amount", "15 μmol", q(15, "umol", "15 μmol"), [ev(4)], "inherited_general_condition")
    set_charge(pcata, "mol_percent", "5 mol%", q(5, "mol%", "5 mol%"), [ev(3), ev(4)])
    pcata["equivalence_reference"] = direct(f"relative to 1{x}", mol(f"1{x}"), e)
    psolv = participant("solvent", mol("dce"), "solvent", e)
    set_charge(psolv, "volume", "3 mL", q(3, "mL", "3 mL"), [ev(3), ev(4)])
    pp = participant("product", mol(f"2{x}"), "product", e, product_role="desired product")
    psilica = participant("silica", mol("silica"), "purification_agent", [ev(4)], stage="purification")
    ppet = participant("petroleum_ether", mol("petether"), "purification_agent", [ev(4)], stage="purification")
    pe = participant("diethyl_ether", mol("et2o"), "purification_agent", [ev(4)], stage="purification")
    r["participants"] = [ps, pcata, psolv, pp, psilica, ppet, pe]
    r["reaction_scale"] = conflict(
        ["Scheme 2: 0.3 mol", "general procedure: 0.3 mmol"],
        [direct("0.3 mol", {"reference_molecule_id": mol(f"1{x}"), "amount": q(0.3, "mol", "0.3 mol")}, [ev(3), ev(3, True)]),
         direct("0.3 mmol", {"reference_molecule_id": mol(f"1{x}"), "amount": q(0.3, "mmol", "0.3 mmol")}, [ev(4)])], e)
    r["solvent_conditions"]["composition"] = direct("DCE", [mol("dce")], e)
    r["solvent_conditions"]["total_volume"] = direct("3 mL", q(3, "mL", "3 mL"), e)
    add_step(r, "reaction", e, ["substrate", "catalyst", "solvent", "product"],
             "gold complex C added to a DCE solution of substrate, then stirred",
             ("room temperature", "room temperature"),
             ("12–24 h", rangeq(12, 24, "h", "12–24 h")))
    add_step(r, "solvent_removal", [ev(4)], ["solvent"], "solvent removed under reduced pressure",
             stage="workup", method="inherited_general_condition")
    add_step(r, "flash_chromatography", [ev(4)], ["silica", "petroleum_ether", "diethyl_ether", "product"],
             "flash column chromatography on silica gel; petroleum ether/Et2O = 100:1 to 10:1",
             stage="purification", method="inherited_general_condition")
    add_outcome(r, f"yield_2{x}", "yield", mol(f"2{x}"), f"{yld}%", q(yld, "%", f"{yld}%"),
                [ev(3), ev(3, True), ev(cp)], role="desired product", yield_type="isolated yield", basis="isolated")
    add_outcome(r, f"amount_2{x}", "product_amount", mol(f"2{x}"), f"{masses[x]} mg", q(masses[x], "mg", f"{masses[x]} mg"),
                [ev(cp), ev(cp, True)], role="desired product", basis="reported product mass")
    if x == "r":
        add_outcome(r, "regioselectivity", "selectivity", mol("2r"), "rr = 10:1", "10:1", [ev(3), ev(6)],
                    scope={"metric": "regioisomeric ratio", "order": "major:minor", "isomer_assignments": "unassigned"})
    if x == "v":
        add_outcome(r, "diastereoselectivity", "selectivity", mol("2v"), "dr = 2.6:1", "2.6:1", [ev(3), ev(6)],
                    scope={"metric": "diastereomeric ratio", "order": "major:minor", "configurations": "unassigned"})
    if x == "w":
        add_outcome(r, "diastereoselectivity", "selectivity", mol("2w"), "dr = 6.1:1", "6.1:1", [ev(3), ev(6)],
                    scope={"metric": "diastereomeric ratio", "order": "major:minor", "configurations": "unassigned"})
    set_classification(r, "C-C", 2, formed=3, evidence=[ev(3), ev(3, True)],
                       note="Sequential cycloisomerization drawing supports three net C-C ring-forming connections; exact atom mapping and accompanying bond-order changes remain unresolved.")
    r["procedure"]["original_text"] = direct(
        "To a solution of 1,3,9-ene-allene-yne 1 (0.3 mmol) in DCE (3 mL) was added gold(I) complex C (14 mg, 15 μmol, 5 mol-%) and the reaction mixture was stirred at room temperature for 12 h. Upon completion, the solvent was removed under reduced pressure and the residue was purified by flash column chromatography on silica gel (eluent: petroleum ether/Et2O = 100:1 to 10:1).",
        "page-4 general procedure", [ev(4)])
    r["procedure"]["general_procedure_reference"] = direct("Experimental Section general procedure for products 2", "page-4 general procedure", [ev(4)])
    r["procedure"]["addition_order"] = direct("substrate in DCE, then gold(I) complex C", [mol(f"1{x}"), mol("dce"), mol("catC")], [ev(4)])
    r["notes"].append(direct("Scheme 2 reports 12–24 h; the page-4 general procedure states 12 h. The exact duration for this substrate is not selected beyond the Scheme 2 range.",
                                     "scope/general-procedure duration specificity unresolved", [ev(3), ev(4)]))
    r["review"]["unresolved"] = [
        "Scheme 2 prints 0.3 mol, whereas the general procedure prints 0.3 mmol; both are preserved as a conflict.",
        "Exact substrate-specific time within 12–24 h is not reported in the provided pages.",
        "Complex reactant/product structures and stereochemical atom assignments require expert review; no concrete SMILES or complete atom mapping is asserted.",
    ]
    reactions.append(r)

# Scheme 3: scale-up and two derivatization arrows.
e3 = [ev(3), ev(3, True)]
r = reaction("scheme3_scaleup_1e_to_2e", e3, categories=["scale_up"], location="Scheme 3, first arrow",
             name="Gram-scale gold-catalyzed sequential cycloisomerization", transformation="sequential cycloisomerization/formal [4+2] cycloaddition")
ps = participant("substrate", mol("1e"), "reactant", e3); set_charge(ps, "amount", "2.1 mmol", q(2.1, "mmol", "2.1 mmol"), e3); set_charge(ps, "mass", "1.04 g", q(1.04, "g", "1.04 g"), e3)
pc = participant("catalyst", mol("catC"), "catalyst", e3); set_charge(pc, "mol_percent", "2 mol%", q(2, "mol%", "2 mol%"), e3); pc["equivalence_reference"] = direct("relative to 1e", mol("1e"), e3)
pv = participant("solvent", mol("dce"), "solvent", e3)
pp = participant("product", mol("2e"), "product", e3, product_role="desired product")
r["participants"] = [ps, pc, pv, pp]
r["reaction_scale"] = direct("1e (2.1 mmol, 1.04 g)", {"reference_molecule_id": mol("1e"), "amount": q(2.1, "mmol"), "mass": q(1.04, "g")}, e3)
r["solvent_conditions"]["composition"] = direct("DCE", [mol("dce")], e3)
add_step(r, "reaction", e3, ["substrate", "catalyst", "solvent", "product"], "stirred", ("rt", "room temperature"), ("12 h", q(12, "h")))
add_outcome(r, "yield_2e", "yield", mol("2e"), "88%", q(88, "%"), e3, role="desired product")
add_outcome(r, "amount_2e", "product_amount", mol("2e"), "0.92 g", q(0.92, "g"), e3, role="desired product")
set_classification(r, "C-C", 2, formed=3, evidence=e3, note="Same sequential carbocyclization class as Scheme 2; exact atom mapping unresolved.")
r["review"]["unresolved"] = ["Yield measurement basis is not explicitly stated for the scale-up arrow.", "Complex structures/atom mapping need expert review."]
reactions.append(r)

r = reaction("scheme3_2e_to_3e", e3, categories=["derivatization"], location="Scheme 3, second arrow",
             name="Diepoxidation of 2e", transformation="double epoxidation")
p1 = participant("substrate", mol("2e"), "reactant", e3)
pr = participant("mcpba", mol("mcpba"), "reagent", e3); set_charge(pr, "equivalents", "3 equiv", q(3, "equiv"), e3); pr["equivalence_reference"] = direct("relative to 2e", mol("2e"), e3)
ps = participant("solvent", mol("dcm"), "solvent", e3)
pp = participant("product", mol("3e"), "product", e3, product_role="diepoxide product")
r["participants"] = [p1, pr, ps, pp]
r["solvent_conditions"]["composition"] = direct("DCM", [mol("dcm")], e3)
add_step(r, "reaction", e3, ["substrate", "mcpba", "solvent", "product"], "treated with m-CPBA", ("rt", "room temperature"), ("24 h", q(24, "h")))
add_outcome(r, "yield_3e", "yield", mol("3e"), "88%", q(88, "%"), e3, role="diepoxide product")
add_outcome(r, "diastereoselectivity", "selectivity", mol("3e"), "single diastereomer", "single diastereomer", e3,
            scope={"metric": "diastereoselectivity", "configuration": "relative drawing shown; absolute assignment unreported"})
set_classification(r, "C-O", 3, formed=4, evidence=e3, note="Two alkene epoxidations visibly create four C-O sigma bonds; atom identities and hydrogen changes unresolved.")
r["review"]["unresolved"] = ["Yield and 'single diastereomer' measurement methods are not reported.", "Absolute configuration is not assigned."]
reactions.append(r)

r = reaction("scheme3_3e_to_4e", e3, categories=["derivatization"], location="Scheme 3, third arrow",
             name="FeCl3-mediated conversion of diepoxide 3e to diketone 4e", transformation="epoxide-to-diketone rearrangement")
p1 = participant("substrate", mol("3e"), "reactant", e3)
pc = participant("fecl3", mol("fecl3"), "catalyst", e3); set_charge(pc, "mol_percent", "10 mol%", q(10, "mol%"), e3); pc["equivalence_reference"] = direct("relative to 3e", mol("3e"), e3)
ps = participant("solvent", mol("dce"), "solvent", e3)
pp = participant("product", mol("4e"), "product", e3, product_role="diketone product")
r["participants"] = [p1, pc, ps, pp]
r["solvent_conditions"]["composition"] = direct("DCE", [mol("dce")], e3)
add_step(r, "reaction", e3, ["substrate", "fecl3", "solvent", "product"], "treated with FeCl3", ("80 °C", q(80, "°C")), ("2 h", q(2, "h")))
add_outcome(r, "yield_4e", "yield", mol("4e"), "68%", q(68, "%"), e3, role="diketone product")
set_classification(r, "DB-C-O", 14, formed=2, change_type="order_increased", before=1, after=2, evidence=e3,
                   note="Two epoxide C-O bonds are promoted to carbonyl C=O bonds in the depicted diketone; accompanying C-O cleavages are not atom-mapped.")
r["review"]["unresolved"] = ["Yield measurement basis is not explicitly stated.", "Complete rearrangement atom mapping and stereochemical fate are unresolved."]
reactions.append(r)

# Scheme 4 controls.
r = reaction("scheme4a_intermediate_not_isolable", e3, categories=["control", "mechanistic_study"], location="Scheme 4a",
             name="Attempted isolation of cyclopentadiene intermediate 5", transformation="1,3-ene-allene cycloisomerization")
p1 = participant("substrate", mol("1y"), "reactant", e3)
pc = participant("catalyst", mol("catC"), "catalyst", e3); set_charge(pc, "mol_percent", "5 mol%", q(5, "mol%"), e3); pc["equivalence_reference"] = direct("relative to 1y", mol("1y"), e3)
ps = participant("solvent", mol("dce"), "solvent", e3)
pp = participant("product", mol("5"), "product", e3, product_role="transient cyclopentadiene intermediate")
r["participants"] = [p1, pc, ps, pp]
r["solvent_conditions"]["composition"] = direct("DCE", [mol("dce")], e3)
add_step(r, "reaction", e3, ["substrate", "catalyst", "solvent", "product"], "gold-catalyzed cycloisomerization", ("rt", "room temperature"), ("1 h", q(1, "h")))
add_step(r, "attempted_chromatography", e3, ["product"], "attempted column chromatography; rapid polymerization observed", stage="purification")
add_outcome(r, "decomposition_5", "decomposition", mol("5"), "5 not isolable; rapid polymerization upon column chromatography", "rapid polymerization; not isolable", e3, role="transient intermediate")
set_classification(r, "C-C", 2, formed=1, evidence=e3, note="Cyclopentadiene ring closure gives at least one new C-C connection; full bond-order rearrangement and atom mapping unresolved.")
r["review"]["unresolved"] = ["No isolated yield or analytical amount for 5 is reported.", "Polymer identity and amount are not reported."]
reactions.append(r)

r = reaction("scheme4b_trapping_to_6", e3, categories=["control", "mechanistic_study"], location="Scheme 4b",
             name="Cyclopentadiene trapping with N-phenylmaleimide", transformation="cycloisomerization followed by intermolecular [4+2] cycloaddition")
p1 = participant("substrate", mol("1y"), "reactant", e3)
pc = participant("catalyst", mol("catC"), "catalyst", e3); set_charge(pc, "mol_percent", "5 mol%", q(5, "mol%"), e3); pc["equivalence_reference"] = direct("relative to 1y", mol("1y"), e3)
ps = participant("solvent", mol("dce"), "solvent", e3)
pr = participant("maleimide", mol("npm"), "reagent", e3); set_charge(pr, "equivalents", "2.0 equiv", q(2.0, "equiv"), e3); pr["equivalence_reference"] = direct("relative to 1y", mol("1y"), e3)
pp = participant("product", mol("6"), "product", e3, product_role="[4+2] trapping adduct")
r["participants"] = [p1, pc, ps, pr, pp]
r["solvent_conditions"]["composition"] = direct("DCE", [mol("dce")], e3)
add_step(r, "cycloisomerization", e3, ["substrate", "catalyst", "solvent"], "gold-catalyzed cycloisomerization", ("rt", "room temperature"), ("1 h", q(1, "h")))
add_step(r, "maleimide_trapping", e3, ["maleimide", "product"], "N-phenylmaleimide added, then reaction continued", ("rt", "room temperature"), ("15 h", q(15, "h")))
add_outcome(r, "yield_6", "yield", mol("6"), "81%", q(81, "%"), e3, role="[4+2] trapping adduct")
set_classification(r, "C-C", 2, formed=3, evidence=e3, note="One cycloisomerization ring closure plus two Diels-Alder C-C connections are inferred from the depicted sequence; atom mapping unresolved.")
r["review"]["unresolved"] = ["Yield measurement basis is not explicitly stated.", "Absolute/relative stereochemical assignments of adduct 6 are not textually reported."]
reactions.append(r)

# Scheme 5 is an author-proposed mechanism, not an experimental execution.
e4 = [ev(4), ev(4, True)]
r = reaction("scheme5_proposed_mechanism", e4, record_type="proposed_mechanism", categories=[], location="Scheme 5",
             name="Proposed mechanism for gold(I)-catalyzed sequential cycloisomerization of 1a",
             transformation="proposed Nazarov cyclization/1,2-hydride shift/formal [4+2] sequence")
r["participants"] = [participant("substrate", mol("1a"), "reactant", e4),
                     participant("gold", mol("generic_au"), "catalyst", e4),
                     participant("product", mol("2a"), "product", e4, product_role="illustrated final product")]
r["notes"] = [direct(
    "Proposed sequence: allene activation (Ia), pentadienyl cation IIa, IIIa/IIIa' gold carbenoid, 1,2-hydride migration/deauration to cyclopentadiene IVa, alkyne activation Va, 5-exo-dig species VIa, formal [4+2] to VIIa, rearomatization/protodeauration to 2a.",
    "author-proposed mechanism; not independent experimental proof", e4)]
set_classification(r, "C-C", 2, formed=3, evidence=e4, note="Overall proposed pathway depicts three net C-C ring-forming connections; elementary bond-order/hydrogen changes remain unresolved.")
r["review"]["unresolved"] = ["Proposed intermediates are not treated as observed compounds.", "No kinetic, isotope, or independently isolated mechanistic evidence is reported on page 4."]
reactions.append(r)

bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:figure1", "kind": "nonreaction_figure",
     "evidence_ids": [ev(1), ev(1, True)],
     "summary": "Figure 1 shows selected natural products as motivation; it is not a synthesis experiment in this paper."},
    {"context_id": f"{DOC}:context:characterization_p4", "kind": "characterization",
     "evidence_ids": [ev(4), ev(4, True)],
     "summary": "Physical descriptions, melting points, NMR and HRMS for 2a-2f are linked context for preparations, not separate reactions."},
    {"context_id": f"{DOC}:context:characterization_p5", "kind": "characterization",
     "evidence_ids": [ev(5), ev(5, True)],
     "summary": "Continuation of 2f and physical/NMR/HRMS data for 2g-2p; these are not separate synthesis records."},
    {"context_id": f"{DOC}:context:characterization_p6", "kind": "characterization",
     "evidence_ids": [ev(6), ev(6, True)],
     "summary": "Continuation of 2p and physical/NMR/HRMS data for 2q-2x; these are not separate synthesis records."},
]

# Persist graph diagnostics and literal mention hints as raw artifacts.
stereo_checks = {}
for m in bundle["molecules"]:
    smi = st.value(m["smiles"])
    if smi:
        stereo_checks[m["molecule_id"]] = st.stereo_inspect(smi)
st.write(RUN / "raw_model/stereo_checks.json", stereo_checks)
st.write(RUN / "raw_model/source_mentions.json", st.source_mentions(DATASET, ["1a", "2a", "1e", "2e", "1y", "3e", "4e", "6"], 180))
st.write(RUN / "submitted.json", bundle)

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "source_context": len(bundle["source_context"]), "stereo_checks": len(stereo_checks)}, ensure_ascii=False))
