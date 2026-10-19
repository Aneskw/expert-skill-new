#!/usr/bin/env python3
"""Source transcription for doc_a92a074d2b5271ac, PDF pages 1-6.

This is a session-local transcription program.  It uses the v3.1 constructors but
does not alter the extraction tools, schemas, or scientific rules.
"""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_a92a074d2b5271ac"
RUN = ROOT / "runs/current"
E = {page: [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"] for page in range(1, 7)}


def f(raw=None, normalized=None, ev=(), method="original_extraction", status=None,
      derivation=None):
    return st.fact(raw, normalized, ev, method=method, status=status,
                   derivation=derivation)


def inferred(normalized, ev, derivation):
    return f(None, normalized, ev, method="inference", derivation=derivation)


def derived(raw, normalized, ev, derivation):
    return f(raw, normalized, ev, method="structure_derivation", derivation=derivation)


def missing(ev=()):
    return f(ev=ev)


def not_applicable():
    return f(method="not_applicable", status="not_applicable")


bundle = json.loads((RUN / "submitted.json").read_text())
molecules = {}


def add_molecule(key, label, name=None, ev=None, smiles=None, stereo="ordinary",
                 note=None, aliases=None):
    ev = ev or E[3]
    mid = f"{DOC}:mol_{key}"
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name or label, name or label, ev)
    m["names"]["english"] = f(name or label, name or label, ev)
    m["names"]["abbreviations"] = f(aliases, aliases, ev) if aliases else missing()
    m["names"]["aliases"] = missing()
    m["cas_numbers"] = missing()
    m["external_properties"] = []
    if smiles:
        diag = st.stereo_inspect(smiles)
        canonical = diag["canonical_isomeric_smiles"]
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            f"structure drawing or named identity for {label}", canonical, ev,
            "Candidate connectivity was transcribed from the supplied source drawing/name and normalized with RDKit; parsing does not establish source correctness.",
        )
        m["generic_template"] = None
        m["stereochemistry"]["isomeric_smiles"] = derived(
            f"candidate structure for {label}", canonical, ev,
            "RDKit canonical isomeric SMILES for the candidate graph; unencoded source stereochemistry remains described separately.",
        )
        m["stereochemistry"]["identity_key"] = derived(
            f"candidate structure for {label}", diag["stereo_graph_key"], ev,
            "SHA-256 stereo graph key calculated by standard stereo inspection; not a sample identity.",
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing()
        m["generic_template"] = None
        m["stereochemistry"]["isomeric_smiles"] = missing()
        m["stereochemistry"]["identity_key"] = missing()
    sc = m["stereochemistry"]
    if stereo == "ordinary":
        sc["atom_centers"] = inferred([], ev, "No tetrahedral stereocenter is present in the candidate graph.") if smiles else missing()
        sc["double_bonds"] = inferred([], ev, "No source-assigned E/Z element applies to the candidate identity.") if smiles else missing()
        sc["axial"] = inferred([], ev, "No configurationally stable axis is present in the candidate graph.") if smiles else missing()
        sc["mixture"] = missing()
        sc["assignment"] = missing()
        sc["chirality_status"] = inferred("achiral", ev, "Candidate graph has no stereogenic center or axis.") if smiles else inferred("unresolved", ev, "The source label is retained but its complete structure was not encoded.")
        sc["geometry_status"] = inferred("not_applicable", ev, "No source-assigned alkene geometry applies.") if smiles else inferred("unresolved", ev, "Complete structure and geometric elements are unresolved.")
        sc["representation_status"] = inferred("not_applicable", ev, "No stereochemical encoding is required for this candidate graph.") if smiles else inferred("unresolved", ev, "No complete structure representation was asserted.")
    elif stereo == "chiral_product":
        sc["atom_centers"] = f("dashed bond at the newly formed allylic stereocenter", [{"source_assignment": "drawn stereocenter; R/S descriptor not printed"}], ev)
        sc["double_bonds"] = inferred([], ev, "The depicted alkene is terminal methylene, so E/Z is not applicable.")
        sc["axial"] = inferred([], ev, "The allene axis is absent from the product.")
        sc["mixture"] = f("enantioenriched product; experiment-specific er is stored on each reaction", "enantioenriched", ev)
        sc["assignment"] = f("configuration drawn; absolute R/S descriptor not printed in supplied pages", {"assignment": "source-drawn", "descriptor": None}, ev)
        sc["chirality_status"] = f("drawn stereocenter", "specified", ev)
        sc["geometry_status"] = inferred("not_applicable", ev, "Terminal methylene has no E/Z descriptor.")
        sc["representation_status"] = inferred("partial" if smiles else "not_encoded", ev, "Connectivity may be represented, but the source-drawn absolute configuration is not safely encoded without the unavailable SI assignment.")
    elif stereo == "two_centers":
        sc["atom_centers"] = f("two drawn stereocenters in Scheme 3", [{"count": 2, "descriptors": None}], ev)
        sc["double_bonds"] = inferred([], ev, "No source-assigned E/Z element applies.")
        sc["axial"] = inferred([], ev, "No stereogenic axis is depicted.")
        sc["mixture"] = f("only diastereomer", "single reported diastereomer", ev)
        sc["assignment"] = f("absolute configuration assigned by X-ray crystallography; R/S descriptors not printed", {"method": "X-ray crystallography", "descriptors": None}, ev)
        sc["chirality_status"] = f("two stereocenters drawn", "specified", ev)
        sc["geometry_status"] = inferred("not_applicable", ev, "Hydrogenated product has no relevant alkene E/Z element.")
        sc["representation_status"] = inferred("partial" if smiles else "not_encoded", ev, "The supplied main text does not print R/S descriptors, so candidate connectivity does not claim complete stereo encoding.")
    elif stereo == "three_centers":
        sc["atom_centers"] = f("Scheme 3 product bearing three stereocenters", [{"count": 3, "descriptors": None}], ev)
        sc["double_bonds"] = missing()
        sc["axial"] = missing()
        sc["mixture"] = f("complete diastereoselectivity / only diastereomer", "single reported diastereomer", ev)
        sc["assignment"] = f("source-drawn relative configuration; absolute descriptors not fully printed", {"descriptors": None}, ev)
        sc["chirality_status"] = f("three stereocenters drawn", "specified", ev)
        sc["geometry_status"] = inferred("not_applicable", ev, "No alkene E/Z assignment is reported for the bicyclic product.")
        sc["representation_status"] = inferred("not_encoded", ev, "Complex bicyclic stereochemistry is retained textually rather than forced into an unsafe SMILES.")
    elif stereo == "linear":
        sc["atom_centers"] = missing()
        sc["double_bonds"] = f("linear allylic regioisomer; E/Z assignment not separately reported", None, ev, status="unreadable")
        sc["axial"] = inferred([], ev, "Allene axis is absent from the product.")
        sc["mixture"] = f("observed as the minor linear regioisomer in the b/l ratio", "regioisomeric mixture component", ev)
        sc["assignment"] = missing()
        sc["chirality_status"] = inferred("unspecified", ev, "No separate chiral assignment is reported for the linear regioisomer.")
        sc["geometry_status"] = inferred("unresolved", ev, "The linear product geometry was not independently assigned in the supplied pages.")
        sc["representation_status"] = inferred("not_encoded", ev, "The minor linear structure is not completely encoded.")
    m["identity_note"] = f(note, note, ev) if note else f(
        "Machine-transcribed source identity; CAS was not reported and structure correctness remains for expert review.",
        "needs expert structure review", ev,
    )
    molecules[key] = m
    return mid


def add_template(key, label, fixed, variables, ev, product=False):
    mid = f"{DOC}:mol_{key}"
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(label, label, ev)
    m["names"]["english"] = f(label, label, ev)
    m["names"]["abbreviations"] = missing()
    m["names"]["aliases"] = missing()
    m["cas_numbers"] = missing()
    m["representation_type"] = "generic_template"
    m["smiles"] = missing()
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, ev),
        "variables": [{
            "label": v[0],
            "attachment_points": f(v[1], v[1], ev),
            "definition": f(v[2], v[2], ev),
            "allowed_values": f(v[3], v[3], ev),
            "constraints": f(v[4], v[4], ev),
        } for v in variables],
        "instance_mapping": [],
        "reason": f("generic R/Ar/n drawing", "The source uses an unexpanded family drawing; no wildcard SMILES is asserted.", ev),
    }
    sc = m["stereochemistry"]
    sc["atom_centers"] = f("drawn product stereocenter" if product else None, [{"descriptor": None}] if product else None, ev) if product else missing()
    sc["double_bonds"] = missing()
    sc["axial"] = missing()
    sc["mixture"] = missing()
    sc["assignment"] = f("enantioenriched product drawing; descriptor not printed", {"descriptor": None}, ev) if product else missing()
    sc["chirality_status"] = inferred("partially_specified" if product else "unresolved", ev, "Generic source drawing does not define every substituent or descriptor.")
    sc["geometry_status"] = inferred("unresolved", ev, "Generic family geometry cannot be fully assigned.")
    sc["representation_status"] = inferred("not_encoded", ev, "Generic R-group structure is retained as a template rather than concrete SMILES.")
    sc["isomeric_smiles"] = missing()
    sc["identity_key"] = missing()
    m["identity_note"] = f("generic source family", "R-group family retained without unsupported expansion", ev)
    m["external_properties"] = []
    molecules[key] = m
    return mid


# Core substrates, common reagents, solvents, and catalyst system.
add_molecule("1a", "1a", "2,3-butadien-1-ol", E[2], "C=C=CCO")
add_molecule("1b", "1b", "3,4-pentadien-1-ol", E[3], "C=C=CCCO")
add_molecule("1c", "1c", "4,5-hexadien-1-ol", E[3], "C=C=CCCCO")
add_molecule("1d", "1d", "n-butyl ether of 2,3-allenol", E[3], "C=C=CCOCCCC")
add_molecule("1e", "1e", "alkyl-substituted allene bearing n-hexyl", E[3], None,
             note="The label and n-hexyl drawing are retained; full connectivity is not forced from the compact scheme.")
aryl_names = {
    "2a": ("iodobenzene", "Ic1ccccc1"),
    "2b": ("4-iodotoluene", "Cc1ccc(I)cc1"),
    "2c": ("3-iodotoluene", "Cc1cccc(I)c1"),
    "2d": ("4-iodoanisole", "COc1ccc(I)cc1"),
    "2e": ("1-fluoro-3-iodobenzene", "Fc1cccc(I)c1"),
    "2f": ("1-bromo-4-iodobenzene", "Brc1ccc(I)cc1"),
    "2g": ("4-iodoacetophenone", "CC(=O)c1ccc(I)cc1"),
    "2h": ("methyl 4-iodobenzoate", "COC(=O)c1ccc(I)cc1"),
    "2i": ("4-iodobenzonitrile", "N#Cc1ccc(I)cc1"),
    "2j": ("1-iodo-4-nitrobenzene", "O=[N+]([O-])c1ccc(I)cc1"),
    "2k": ("2-iodonaphthalene", "Ic1ccc2ccccc2c1"),
    "2l": ("1-iodonaphthalene", "Ic1cccc2ccccc12"),
}
for key, (name, smi) in aryl_names.items():
    add_molecule(key, key, name, E[3], smi)
add_molecule("3a", "3a", "2-pyridone", E[2], "O=c1cccc[nH]1")
pyr_names = {
    "3b": "5-methyl-2-pyridone", "3c": "4-methyl-2-pyridone",
    "3d": "6-methyl-2-pyridone", "3e": "4-chloro-2-pyridone",
    "3f": "5-nitro-2-pyridone", "3g": "5-formyl-2-pyridone",
    "3h": "4-cyano-2-pyridone", "3i": "3,5-dibromo-2-pyridone",
}
for key, name in pyr_names.items():
    add_molecule(key, key, name, E[3], None,
                 note="Substituent identity/position is transcribed from Figure 3; a complete SMILES is withheld pending expert drawing review.")
add_molecule("pdcat", "Pd2(dba)3·CHCl3", "tris(dibenzylideneacetone)dipalladium(0) chloroform adduct", E[2])
for i in range(1, 9):
    add_molecule(f"L{i}", f"L{i}", "(S,S,S)-SKP" if i == 6 else f"chiral ligand L{i}", E[2],
                 note="Ligand structure is source-drawn in Table 1; no complete SMILES was asserted.")
for key, label, name, smi in [
    ("K2CO3", "K2CO3", "potassium carbonate", "[K+].[K+].[O-]C([O-])=O"),
    ("K3PO4", "K3PO4", "potassium phosphate", "[K+].[K+].[K+].[O-]P([O-])([O-])=O"),
    ("Cs2CO3", "Cs2CO3", "cesium carbonate", "[Cs+].[Cs+].[O-]C([O-])=O"),
    ("KOH", "KOH", "potassium hydroxide", "[K+].[OH-]"),
    ("KF", "KF", "potassium fluoride", "[K+].[F-]"),
    ("Et3N", "Et3N", "triethylamine", "CCN(CC)CC"),
    ("toluene", "toluene", "toluene", "Cc1ccccc1"),
    ("THF", "THF", "tetrahydrofuran", "C1CCOC1"),
    ("CPME", "CPME", "cyclopentyl methyl ether", "COC1CCCC1"),
    ("CH2Cl2", "CH2Cl2", "dichloromethane", "ClCCl"),
    ("CH3CN", "CH3CN", "acetonitrile", "CC#N"),
    ("MTBE", "MTBE", "methyl tert-butyl ether", "COC(C)(C)C"),
    ("MeOH", "MeOH", "methanol", "CO"),
    ("H2", "H2", "hydrogen", "[H][H]"),
    ("PdC", "Pd/C", "palladium on carbon", None),
    ("BF3OEt2", "BF3·OEt2", "boron trifluoride diethyl etherate", None),
]:
    add_molecule(key, label, name, E[2] if key not in {"MeOH", "H2", "PdC", "BF3OEt2"} else E[4], smi)
add_molecule("8", "8", "3-nitrobenzaldehyde", E[4], "O=Cc1cccc([N+](=O)[O-])c1")

# Generic R-group families.
gev = E[1] + E[3]
add_template("allenol_family", "allenol 1", "terminal allene-(CH2)n-OH", [
    ("n", "chain length between allene and OH", "integer chain length", [1, 2], "n = 1 or 2 in the successful family"),
], gev)
add_template("aryl_iodide_family", "aryl iodide 2", "Ar-I", [
    ("Ar", "carbon bonded to iodine", "aryl group", "substituted phenyl or naphthyl shown in scope", "must match a source-drawn aryl iodide"),
], gev)
add_template("pyridone_family", "2-pyridone 3", "2-pyridone N-H core", [
    ("R", "pyridone-ring substituent positions", "ring substituent(s)", "H, methyl, nitro, formyl, cyano, halogen", "positions follow source drawings"),
], gev)
add_template("product_family", "branched N-allylic 2-pyridone 4/6", "branched allylic alcohol bearing Ar and N-pyridone substituents", [
    ("Ar", "vinylic aryl position", "aryl group", "scope aryl groups", "inherits the aryl iodide identity"),
    ("R", "pyridone-ring positions", "pyridone substituent(s)", "scope substituents", "inherits the pyridone identity"),
    ("n", "alcohol tether length", "integer chain length", [1, 2], "source family n = 1 or 2"),
], gev, product=True)

# Product identities. 4aa and 7 have connectivity candidates; other complex products remain source-labelled.
add_molecule("4aa", "4aa", "branched aminoarylation product 4aa", E[2] + E[3] + E[4],
             "C=C(c1ccccc1)C(CO)n1ccccc1=O", "chiral_product",
             note="Table 1 states that absolute configuration was assigned by X-ray crystallography, but the supplied main pages do not print an R/S descriptor; candidate connectivity therefore omits an @/@@ claim.")
add_molecule("5aa", "5aa", "linear regioisomer 5aa", E[2], None, "linear")

fig2_products = {
    "4ba": (67, "97.5:2.5", "75:25", "2b"),
    "4ca": (58, "95:5", "72:28", "2c"),
    "4da": (65, "97:3", "75:25", "2d"),
    "4ea": (71, "94:6", "77:23", "2e"),
    "4fa": (75, "97:3", "82:18", "2f"),
    "4ga": (75, "98:2", "81:19", "2g"),
    "4ha": (72, "98:2", "80:20", "2h"),
    "4ia": (65, "98:2", "75:25", "2i"),
    "4ja": (67, "92:8", "75:25", "2j"),
    "4ka": (76, "98:2", "82:18", "2k"),
    "4la": (35, "56:44", "70:30", "2l"),
}
for code in fig2_products:
    add_molecule(code, code, f"branched aminoarylation product {code}", E[3], None, "chiral_product")
    add_molecule("5" + code[1:], "linear regioisomer corresponding to " + code,
                 f"linear regioisomer corresponding to {code}", E[3], None, "linear")

fig3_products = {
    "4ab": (73, "98:2", "82:18", "2a", "3b", "A"),
    "4ac": (58, "92:8", "80:20", "2a", "3c", "A"),
    "4ad": (60, "91:9", "71:29", "2a", "3d", "B"),
    "4ae": (60, "91:9", "70:30", "2a", "3e", "B"),
    "4af": (67, "93:7", "74:26", "2a", "3f", "B"),
    "4ag": (64, "95:5", "72:28", "2a", "3g", "B"),
    "4ah": (58, "91:9", "67:33", "2a", "3h", "B"),
    "4ai": (61, "92:8", "71:29", "2a", "3i", "B"),
    "4fb": (73, "91:9", "80:20", "2f", "3b", "B"),
    "4ff": (62, "97:3", "71:29", "2f", "3f", "B"),
    "4kb": (68, "92:8", "76:24", "2k", "3b", "B"),
    "4kf": (65, "97:3", "72:28", "2k", "3f", "B"),
}
for code in fig3_products:
    add_molecule(code, code, f"branched aminoarylation product {code}", E[3], None, "chiral_product")
    add_molecule("5" + code[1:], "linear regioisomer corresponding to " + code,
                 f"linear regioisomer corresponding to {code}", E[3], None, "linear")

for code in ["6a", "6b", "6c", "6d", "6e"]:
    add_molecule(code, code, f"branched aminoarylation product {code}", E[3], None, "chiral_product")
    add_molecule("lin_" + code, "linear regioisomer corresponding to " + code,
                 f"linear regioisomer corresponding to {code}", E[3], None, "linear")
add_molecule("7", "7", "selectively hydrogenated derivative 7", E[4],
             "CC(c1ccccc1)C(CO)n1ccccc1=O", "two_centers")
add_molecule("9", "9", "bicyclic pyridinium tetrafluoroborate 9", E[4], None, "three_centers")
add_molecule("10", "10", "hydrolyzed chiral 2-pyridone 10", E[4], None, "three_centers")


def charge_part(pid, molkey, role, ev, amount=None, equiv=None, molpct=None,
                volume=None, reference=None, stage="reaction"):
    p = st.new_participant(pid, molecules[molkey]["molecule_id"], role, ev, stage)
    p["product_role"] = missing()
    if amount is not None:
        p["charges"]["amount"] = f(amount, amount, ev)
    if equiv is not None:
        p["charges"]["equivalents"] = f(equiv, equiv, ev)
    if molpct is not None:
        p["charges"]["mol_percent"] = f(molpct, molpct, ev)
    if volume is not None:
        p["charges"]["volume"] = f(volume, volume, ev)
    if reference is not None:
        p["equivalence_reference"] = f(reference, molecules[reference]["molecule_id"], ev)
    return p


def product_part(pid, molkey, ev, role="product"):
    p = st.new_participant(pid, molecules[molkey]["molecule_id"], role, ev)
    p["product_role"] = f("major branched product" if role == "product" and not molkey.startswith(("5", "lin_")) else "minor linear regioisomer",
                          "major" if role == "product" and not molkey.startswith(("5", "lin_")) else "minor", ev)
    return p


def classification(ev, pairs=None, coverage="partial", limitations=None):
    pairs = pairs or []
    changes = []
    labels = {}
    for idx, (a, b, kind) in enumerate(pairs, 1):
        before, after = (0, 1) if kind in {"formed", "hydrogen_gain"} else (1, 0)
        change = {"change_id": f"b{idx}", "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": before, "after_order": after, "change_type": kind}
        if kind.startswith("hydrogen_"):
            change["hydrogen_count"] = 1
        changes.append(change)
        if kind in {"formed", "hydrogen_gain"}:
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(change["change_id"])
    deriv = "Partial net bond-change inventory inferred from the source transformation drawing; no complete atom map or coproduct balance is asserted."
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": missing(),
        "bond_changes": [derived(None, c, ev, deriv) for c in changes],
        "labels": derived(None, list(labels.values()), ev, deriv) if coverage not in {"unresolved", "not_applicable"} else inferred([], ev, "No bond-formation label is asserted for an unresolved or no-reaction transformation."),
        "coverage": inferred(coverage, ev, deriv if coverage == "partial" else "No observed transformation is available for bond classification."),
        "limitations": limitations or ["No complete atom mapping.", "Hydrogen and coproduct balance is unresolved.", "Classification is inferred, not a source claim."],
        "model_votes": [],
    }


def set_common_reaction(r, ev, categories, location, general=None, name="Pd-catalyzed regio- and enantioselective aminoarylation of allenols with aryl iodides and 2-pyridones"):
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = f(categories, categories, ev)
    r["location"] = f(location, location, ev)
    r["general_scheme_ids"] = general or []
    r["reaction_name"] = f(name, name, ev)
    r["transformation_type"] = inferred("three-component allene aminoarylation", ev, "Normalized from the authors' reaction name.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def reaction_step(ev, pids, temp=None, time=None, operation=None):
    s = st.new_step("reaction", ev)
    s["participant_ids"] = pids
    s["operation"] = f(operation, operation, ev) if operation else missing()
    s["temperature"] = f(temp, temp, ev) if temp is not None else missing()
    s["time"] = f(time, time, ev) if time is not None else missing()
    s["atmosphere"] = missing()
    s["pressure"] = missing()
    return s


def outcome(oid, kind, targets, ev, value=None, role=None, yield_type=None,
            basis=None, scope=None):
    mids = [molecules[k]["molecule_id"] for k in targets]
    o = st.new_outcome(oid, kind, mids, ev)
    o["product_role"] = f(role, role, ev) if role is not None else missing()
    o["value"] = f(value, value, ev) if value is not None else missing()
    o["yield_type"] = f(yield_type, yield_type, ev) if yield_type is not None else missing()
    o["measurement_basis"] = f(basis, basis, ev) if basis is not None else missing()
    o["scope"] = f(scope, scope, ev) if scope is not None else missing()
    return o


reactions = []

# One family-level scheme from the abstract graphic and Figure 2/3 headers.
rid_generic = f"{DOC}:scheme_general_aminoarylation"
r = st.new_reaction(rid_generic, DOC, gev)
r["record_type"] = "generic_scheme"
r["experiment_categories"] = f("general reaction family", [], gev)
r["location"] = f("page 1 abstract graphic; page 3 Figures 2 and 3 headers", {"pages": [1, 3]}, gev)
r["reaction_name"] = f("Pd-catalyzed regio- and enantioselective aminoarylation of allenols with aryl iodides and 2-pyridones", "Pd-catalyzed regio- and enantioselective aminoarylation", gev)
r["transformation_type"] = inferred("three-component allene aminoarylation", gev, "Normalized from source title and scheme.")
r["participants"] = [
    charge_part("allenol_family", "allenol_family", "reactant", gev),
    charge_part("aryl_iodide_family", "aryl_iodide_family", "reactant", gev),
    charge_part("pyridone_family", "pyridone_family", "reactant", gev),
    product_part("product_family", "product_family", gev),
]
r["steps"] = [reaction_step(gev, [p["participant_id"] for p in r["participants"]], operation="three-component aminoarylation")]
r["classification"] = classification(gev, [("C", "C", "formed"), ("C", "N", "formed")])
r["evidence_ids"] = gev
r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
               "unresolved": ["Generic R/Ar/n families are not complete structures or atom maps."]}
r["completeness"] = "partial"
reactions.append(r)

# Scheme 2 has a distinct generic arrow and includes both successful and failed allene variants.
rid_scheme2 = f"{DOC}:scheme2_allene_variants"
r = st.new_reaction(rid_scheme2, DOC, E[3])
r["record_type"] = "generic_scheme"
r["experiment_categories"] = f("Scheme 2 general arrow", [], E[3])
r["location"] = f("page 3, Scheme 2 general arrow", {"page": 3, "region": "Scheme 2"}, E[3])
r["reaction_name"] = f("Pd-catalyzed aminoarylation of allenes 1 with 2 and 3", "Pd-catalyzed aminoarylation", E[3])
r["transformation_type"] = inferred("three-component allene aminoarylation", E[3], "Normalized from Scheme 2 title.")
r["participants"] = [charge_part("allene1", "allenol_family", "reactant", E[3]), charge_part("aryl2", "aryl_iodide_family", "reactant", E[3]), charge_part("pyridone3", "pyridone_family", "reactant", E[3]), product_part("product6", "product_family", E[3])]
r["steps"] = [reaction_step(E[3], [p["participant_id"] for p in r["participants"]], {"number": 60, "unit": "°C"}, None, "aminoarylation under Scheme 2 conditions")]
r["classification"] = classification(E[3], [("C", "C", "formed"), ("C", "N", "formed")])
r["evidence_ids"] = E[3]
r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [], "unresolved": ["Scheme 2 is a generic R-group drawing; concrete rows are recorded separately."]}
r["completeness"] = "partial"
reactions.append(r)

# Proposed Heck insertion / allylpalladium / hydrogen-bonded substitution sequence on page 1.
r = st.new_reaction(f"{DOC}:proposed_mechanism", DOC, E[1])
r["record_type"] = "proposed_mechanism"
r["experiment_categories"] = f("author-proposed mechanistic sequence", [], E[1])
r["location"] = f("page 1, Scheme 1c, intermediates I and II", {"page": 1, "region": "Scheme 1c"}, E[1])
r["reaction_name"] = f("Heck-type arylation followed by allylic substitution", "proposed aminoarylation mechanism", E[1])
r["transformation_type"] = inferred("proposed organopalladium mechanism", E[1], "Author discussion proposes intermediates I and II; this is not an isolated-intermediate experiment.")
r["participants"] = [charge_part("allenol", "allenol_family", "reactant", E[1]), charge_part("aryl", "aryl_iodide_family", "reactant", E[1]), charge_part("pyridone", "pyridone_family", "reactant", E[1]), product_part("product", "product_family", E[1])]
r["classification"] = classification(E[1], [("C", "C", "formed"), ("C", "N", "formed")])
r["notes"] = [f("authors propose allylpalladium intermediate I and hydrogen-bonded intermediate II", "proposed, not experimentally isolated", E[1])]
r["evidence_ids"] = E[1]
r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [], "unresolved": ["Proposed intermediates are not experimental products."]}
r["completeness"] = "partial"
reactions.append(r)

# Table 1: every row is an independent optimization experiment.
table_rows = {
    1: ("L1", "K2CO3", 2, "toluene", 10, "85:15", "71:29"),
    2: ("L2", "K2CO3", 2, "toluene", 48, "73:27", "75:25"),
    3: ("L3", "K2CO3", 2, "toluene", 30, "90:10", "74:26"),
    4: ("L4", "K2CO3", 2, "toluene", None, None, None),
    5: ("L5", "K2CO3", 2, "toluene", 50, "77:23", "67:33"),
    6: ("L6", "K2CO3", 2, "toluene", 73, "93:7", "80:20"),
    7: ("L7", "K2CO3", 2, "toluene", None, None, None),
    8: ("L8", "K2CO3", 2, "toluene", None, None, None),
    9: ("L6", "K2CO3", 1, "toluene", 70, "98:2", "79:21"),
    10: ("L6", "K2CO3", 1, "THF", 59, "81:19", "67:33"),
    11: ("L6", "K2CO3", 1, "CPME", 67, "92:8", "74:26"),
    12: ("L6", "K2CO3", 1, "CH2Cl2", 45, "70:30", "57:43"),
    13: ("L6", "K2CO3", 1, "CH3CN", 60, "88:12", "71:29"),
    14: ("L6", "Cs2CO3", 1, "toluene", 68, "96:4", "78:22"),
    15: ("L6", "K3PO4", 1, "toluene", 60, "97:3", "74:26"),
    16: ("L6", "KOH", 1, "toluene", 73, "73:27", "76:24"),
    17: ("L6", "KF", 1, "toluene", 30, "93:7", "60:40"),
    18: ("L6", "Et3N", 1, "toluene", 42, "93:7", "60:40"),
}
for entry, (lig, base, eq, solv, yld, er, bl) in table_rows.items():
    ev = E[2] + (E[3] if entry == 9 else [])
    rid = f"{DOC}:table1_entry_{entry:02d}"
    r = set_common_reaction(st.new_reaction(rid, DOC, ev), ev,
                            ["optimization"] + (["substrate_scope"] if entry == 9 else []),
                            {"page": 2, "table": 1, "entry": entry}, [rid_generic])
    base_amount = {"number": 0.4 if eq == 2 else 0.2, "unit": "mmol"}
    r["participants"] = [
        charge_part("p_1a", "1a", "reactant", ev, amount={"number": 0.2, "unit": "mmol"}),
        charge_part("p_2a", "2a", "reactant", ev, amount={"number": 0.24, "unit": "mmol"}),
        charge_part("p_3a", "3a", "reactant", ev, amount={"number": 0.24, "unit": "mmol"}),
        charge_part("p_pd", "pdcat", "catalyst", ev, molpct={"number": 2.5, "unit": "mol%"}, reference="1a"),
        charge_part("p_lig", lig, "ligand", ev, molpct={"number": 5, "unit": "mol%"}, reference="1a"),
        charge_part("p_base", base, "base" if False else "additive", ev, amount=base_amount, equiv={"number": eq, "unit": "equiv"}, reference="1a"),
        charge_part("p_solvent", solv, "solvent", ev, volume={"number": 2.0, "unit": "mL"}),
    ]
    if yld is not None:
        r["participants"] += [product_part("p_4aa", "4aa", ev), product_part("p_5aa", "5aa", ev)]
    r["reaction_scale"] = f("1a (0.2 mmol)", {"reference_molecule_id": molecules["1a"]["molecule_id"], "amount": {"number": 0.2, "unit": "mmol"}}, ev)
    r["solvent_conditions"]["composition"] = f(solv, [molecules[solv]["molecule_id"]], ev)
    r["solvent_conditions"]["total_volume"] = f("2.0 mL", {"number": 2.0, "unit": "mL"}, ev)
    pids = [p["participant_id"] for p in r["participants"]]
    r["steps"] = [reaction_step(ev, pids, {"number": 60, "unit": "°C"}, {"number": 16, "unit": "h"}, "reaction under Table 1 conditions")]
    r["procedure"]["original_text"] = f("Pd2(dba)3·CHCl3 (2.5 mol%), ligand (5 mol%), base (0.2 or 0.4 mmol), 1a (0.2 mmol), 2a (0.24 mmol), 3a (0.24 mmol), solvent (2.0 mL), 60 °C, 16 h", "Table 1 footnote a", ev)
    if yld is None:
        r["outcomes"] = [outcome("nr", "no_reaction", [], ev, "NR", scope={"table_column": "yield"})]
        r["classification"] = classification(ev, coverage="not_applicable")
        r["review"]["unresolved"] = ["Table reports NR; analytical method and material balance are not stated."]
    else:
        r["outcomes"] = [
            outcome("yield", "yield", ["4aa"], ev, {"number": yld, "unit": "%"}, "branched product 4aa", "isolated", "isolated yield", {"table": 1, "entry": entry}),
            outcome("er", "selectivity", ["4aa"], ev, er + " er", "branched product 4aa", None, "HPLC on a chiral stationary phase", {"metric": "er", "order": "major:minor"}),
            outcome("bl", "selectivity", ["4aa", "5aa"], ev, bl + " 4aa:5aa", "branched/linear products", None, "1H NMR spectroscopy of the crude reaction mixture", {"metric": "branched/linear", "order": "4aa:5aa"}),
        ]
        r["classification"] = classification(ev, [("C", "C", "formed"), ("C", "N", "formed")])
        r["review"]["unresolved"] = ["Workup, purification details, atmosphere, and pressure are unavailable without SI.", "Absolute R/S descriptor for 4aa is not printed on supplied pages."]
    r["evidence_ids"] = ev
    reactions.append(r)


def scope_reaction(rid_suffix, code, aryl, pyr, yld, er, bl, condition, ev, general_id=rid_generic):
    rid = f"{DOC}:{rid_suffix}"
    r = set_common_reaction(st.new_reaction(rid, DOC, ev), ev, ["substrate_scope"],
                            {"page": 3, "figure_or_scheme": rid_suffix}, [general_id])
    if condition == "A":
        base, solv, temp, tm = "K2CO3", "toluene", 60, 16
    else:
        base, solv, temp, tm = "K3PO4", "MTBE", 30, 24
    lin = "5" + code[1:]
    r["participants"] = [
        charge_part("p_allenol", "1a", "reactant", ev, amount={"number": 0.2, "unit": "mmol"}),
        charge_part("p_aryl", aryl, "reactant", ev, amount={"number": 0.24, "unit": "mmol"}),
        charge_part("p_pyr", pyr, "reactant", ev, amount={"number": 0.24, "unit": "mmol"}),
        charge_part("p_pd", "pdcat", "catalyst", ev, molpct={"number": 2.5, "unit": "mol%"}, reference="1a"),
        charge_part("p_lig", "L6", "ligand", ev, molpct={"number": 5, "unit": "mol%"}, reference="1a"),
        charge_part("p_base", base, "additive", ev, amount={"number": 0.2, "unit": "mmol"}, equiv={"number": 1, "unit": "equiv"}, reference="1a"),
        charge_part("p_solvent", solv, "solvent", ev, volume={"number": 2.0, "unit": "mL"}),
        product_part("p_major", code, ev), product_part("p_linear", lin, ev),
    ]
    r["reaction_scale"] = f("1a (0.2 mmol)", {"reference_molecule_id": molecules["1a"]["molecule_id"], "amount": {"number": 0.2, "unit": "mmol"}}, ev)
    r["solvent_conditions"]["composition"] = f(solv, [molecules[solv]["molecule_id"]], ev)
    r["solvent_conditions"]["total_volume"] = f("2.0 mL", {"number": 2.0, "unit": "mL"}, ev)
    r["steps"] = [reaction_step(ev, [p["participant_id"] for p in r["participants"]], {"number": temp, "unit": "°C"}, {"number": tm, "unit": "h"}, f"reaction under condition {condition}")]
    r["procedure"]["general_procedure_reference"] = f(f"Figure 2 conditions / condition {condition}", condition, ev)
    r["outcomes"] = [
        outcome("yield", "yield", [code], ev, {"number": yld, "unit": "%"}, "branched product", "isolated", "isolated material", {"product": code}),
        outcome("er", "selectivity", [code], ev, er + " er", "branched product", None, "HPLC on a chiral stationary phase", {"metric": "er", "order": "major:minor"}),
        outcome("bl", "selectivity", [code, lin], ev, bl + " b:l", "branched/linear products", None, "1H NMR spectroscopy of the crude reaction mixture", {"metric": "branched/linear", "order": "branched:linear"}),
    ]
    r["classification"] = classification(ev, [("C", "C", "formed"), ("C", "N", "formed")])
    r["evidence_ids"] = ev
    r["review"]["unresolved"] = ["Workup, purification, atmosphere, and pressure are unreported on supplied pages.", "Product absolute descriptor is not printed; source-drawn stereochemistry is not encoded."]
    return r


# Figure 2: 4aa is already Table 1 entry 9; add the eleven nonduplicate scope rows.
for code, (yld, er, bl, aryl) in fig2_products.items():
    reactions.append(scope_reaction("figure2_" + code, code, aryl, "3a", yld, er, bl, "A", E[3]))

# Figure 3: twelve independent scope rows.
for code, (yld, er, bl, aryl, pyr, cond) in fig3_products.items():
    reactions.append(scope_reaction("figure3_" + code, code, aryl, pyr, yld, er, bl, cond, E[3]))

# Scheme 2: one success for 1b, three NR controls, four additional 1b scope successes.
scheme2_success = {
    "6a": ("1b", "2a", "3a", 81, "94:6", "86:14"),
    "6b": ("1b", "2a", "3b", 82, "90:10", "89:11"),
    "6c": ("1b", "2a", "3g", 84, "92:8", "89:11"),
    "6d": ("1b", "2d", "3a", 78, "91:9", "86:14"),
    "6e": ("1b", "2f", "3a", 86, "89:11", "90:10"),
}
for code, (allene, aryl, pyr, yld, er, bl) in scheme2_success.items():
    ev = E[3]
    rid = f"{DOC}:scheme2_{code}"
    r = set_common_reaction(st.new_reaction(rid, DOC, ev), ev, ["substrate_scope"],
                            {"page": 3, "scheme": 2, "product": code}, [rid_scheme2])
    lin = "lin_" + code
    r["participants"] = [
        charge_part("p_allene", allene, "reactant", ev, amount={"number": 0.2, "unit": "mmol"}),
        charge_part("p_aryl", aryl, "reactant", ev, amount={"number": 0.24, "unit": "mmol"}),
        charge_part("p_pyr", pyr, "reactant", ev, amount={"number": 0.24, "unit": "mmol"}),
        charge_part("p_pd", "pdcat", "catalyst", ev, molpct={"number": 2.5, "unit": "mol%"}, reference=allene),
        charge_part("p_lig", "L6", "ligand", ev, molpct={"number": 5, "unit": "mol%"}, reference=allene),
        charge_part("p_base", "K2CO3", "additive", ev, amount={"number": 0.2, "unit": "mmol"}, equiv={"number": 1, "unit": "equiv"}, reference=allene),
        charge_part("p_solvent", "toluene", "solvent", ev, volume={"number": 2.0, "unit": "mL"}),
        product_part("p_major", code, ev), product_part("p_linear", lin, ev),
    ]
    r["reaction_scale"] = f("optimized conditions A (Table 1 entry 9)", {"reference_molecule_id": molecules[allene]["molecule_id"], "amount": {"number": 0.2, "unit": "mmol"}}, ev, method="inherited_general_condition")
    r["solvent_conditions"]["composition"] = f("toluene", [molecules["toluene"]["molecule_id"]], ev, method="inherited_general_condition")
    r["solvent_conditions"]["total_volume"] = f("2.0 mL under condition A", {"number": 2.0, "unit": "mL"}, ev, method="inherited_general_condition")
    r["steps"] = [reaction_step(ev, [p["participant_id"] for p in r["participants"]], {"number": 60, "unit": "°C"}, {"number": 16, "unit": "h"}, "optimized conditions A")]
    r["procedure"]["general_procedure_reference"] = f("under optimized conditions A (Table 1, entry 9)", "condition A", ev)
    r["outcomes"] = [
        outcome("yield", "yield", [code], ev, {"number": yld, "unit": "%"}, "branched product", None, None, {"product": code}),
        outcome("er", "selectivity", [code], ev, er + " er", "branched product", None, None, {"metric": "er", "order": "major:minor"}),
        outcome("bl", "selectivity", [code, lin], ev, bl + " b:l", "branched/linear products", None, None, {"metric": "branched/linear", "order": "branched:linear"}),
    ]
    r["classification"] = classification(ev, [("C", "C", "formed"), ("C", "N", "formed")])
    r["evidence_ids"] = ev
    r["review"]["unresolved"] = ["Scheme 2 does not state yield, er, or b/l measurement methods.", "Condition A quantities/time are inherited from the explicit cross-reference; workup and purification remain unreported."]
    reactions.append(r)

for substrate in ["1c", "1d", "1e"]:
    ev = E[3]
    rid = f"{DOC}:scheme2_{substrate}_NR"
    r = set_common_reaction(st.new_reaction(rid, DOC, ev), ev, ["substrate_scope", "control"],
                            {"page": 3, "scheme": 2, "substrate": substrate}, [rid_scheme2])
    r["participants"] = [
        charge_part("p_allene", substrate, "reactant", ev, amount={"number": 0.2, "unit": "mmol"}),
        charge_part("p_aryl", "2a", "reactant", ev, amount={"number": 0.24, "unit": "mmol"}),
        charge_part("p_pyr", "3a", "reactant", ev, amount={"number": 0.24, "unit": "mmol"}),
        charge_part("p_pd", "pdcat", "catalyst", ev, molpct={"number": 2.5, "unit": "mol%"}, reference=substrate),
        charge_part("p_lig", "L6", "ligand", ev, molpct={"number": 5, "unit": "mol%"}, reference=substrate),
        charge_part("p_base", "K2CO3", "additive", ev, amount={"number": 0.2, "unit": "mmol"}, equiv={"number": 1, "unit": "equiv"}, reference=substrate),
        charge_part("p_solvent", "toluene", "solvent", ev, volume={"number": 2.0, "unit": "mL"}),
    ]
    r["reaction_scale"] = f("optimized conditions A", {"reference_molecule_id": molecules[substrate]["molecule_id"], "amount": {"number": 0.2, "unit": "mmol"}}, ev, method="inherited_general_condition")
    r["solvent_conditions"]["composition"] = f("toluene", [molecules["toluene"]["molecule_id"]], ev, method="inherited_general_condition")
    r["solvent_conditions"]["total_volume"] = f("2.0 mL", {"number": 2.0, "unit": "mL"}, ev, method="inherited_general_condition")
    r["steps"] = [reaction_step(ev, [p["participant_id"] for p in r["participants"]], {"number": 60, "unit": "°C"}, {"number": 16, "unit": "h"}, "optimized conditions A")]
    r["outcomes"] = [outcome("nr", "no_reaction", [], ev, "NR", scope={"substrate": substrate})]
    r["classification"] = classification(ev, coverage="not_applicable")
    r["evidence_ids"] = ev
    r["review"]["unresolved"] = ["NR observation method, recovery, and material balance are unreported."]
    reactions.append(r)

# Scheme 3a: gram-scale repetition. Preserve the printed "As Table 2" reference without silently correcting it.
ev = E[4]
r = set_common_reaction(st.new_reaction(f"{DOC}:scheme3a_scale_up", DOC, ev), ev, ["scale_up"], {"page": 4, "scheme": "3a"}, [rid_generic])
r["participants"] = [
    charge_part("p_1a", "1a", "reactant", ev, amount={"number": 7, "unit": "mmol"}),
    charge_part("p_2a", "2a", "reactant", ev, amount={"number": 8.4, "unit": "mmol"}),
    charge_part("p_3a", "3a", "reactant", ev, amount={"number": 8.4, "unit": "mmol"}),
    product_part("p_4aa", "4aa", ev), product_part("p_5aa", "5aa", ev),
]
r["reaction_scale"] = f("1a (7 mmol)", {"reference_molecule_id": molecules["1a"]["molecule_id"], "amount": {"number": 7, "unit": "mmol"}}, ev)
r["steps"] = [reaction_step(ev, [p["participant_id"] for p in r["participants"]], operation="three-component aminoarylation under standard conditions")]
r["procedure"]["general_procedure_reference"] = f("As Table 2", "As Table 2", ev)
r["outcomes"] = [
    outcome("yield", "yield", ["4aa"], ev, {"number": 71, "unit": "%"}, "branched product 4aa", None, None, {"scale": "7 mmol 1a"}),
    outcome("mass", "product_amount", ["4aa"], ev, {"number": 1.2, "unit": "g"}, "branched product 4aa", None, None, {"scale": "7 mmol 1a"}),
    outcome("er", "selectivity", ["4aa"], ev, "96:4 er", "branched product 4aa", None, None, {"metric": "er", "order": "major:minor"}),
    outcome("bl", "selectivity", ["4aa", "5aa"], ev, "79:21 b:l", "branched/linear products", None, None, {"metric": "branched/linear", "order": "4aa:5aa"}),
]
r["classification"] = classification(ev, [("C", "C", "formed"), ("C", "N", "formed")])
r["evidence_ids"] = ev
r["review"]["unresolved"] = ["Scheme prints 'As Table 2', while prose says 'standard conditions'; no Table 2 is supplied and the reference is not silently corrected.", "Catalyst/base/solvent/temperature/time, yield basis, and selectivity methods are not restated for this experiment."]
r["notes"] = [f(["As Table 2", "under standard conditions"], None, ev, status="conflict")]
r["notes"][0]["alternatives"] = [f("As Table 2", "As Table 2", ev), f("under standard conditions", "standard conditions", ev)]
reactions.append(r)

# Scheme 3b: selective alkene hydrogenation.
r = set_common_reaction(st.new_reaction(f"{DOC}:scheme3b_hydrogenation", DOC, ev), ev, ["derivatization"], {"page": 4, "scheme": "3b"}, [], "selective hydrogenation of 4aa")
r["transformation_type"] = inferred("alkene hydrogenation", ev, "Normalized from the author's 'selective hydrogenation' description.")
r["participants"] = [charge_part("p_4aa", "4aa", "reactant", ev), charge_part("p_PdC", "PdC", "catalyst", ev), charge_part("p_H2", "H2", "reagent", ev), charge_part("p_MeOH", "MeOH", "solvent", ev), product_part("p_7", "7", ev)]
r["solvent_conditions"]["composition"] = f("MeOH", [molecules["MeOH"]["molecule_id"]], ev)
r["steps"] = [reaction_step(ev, [p["participant_id"] for p in r["participants"]], operation="selective hydrogenation")]
r["outcomes"] = [outcome("yield", "yield", ["7"], ev, {"number": 95, "unit": "%"}, "product 7", None, None, {"scheme": "3b"}), outcome("dr", "selectivity", ["7"], ev, "only diastereomer", "product 7", None, None, {"metric": "diastereoselectivity"})]
r["classification"] = classification(ev, [("C", "H", "hydrogen_gain"), ("C", "H", "hydrogen_gain")], limitations=["Net addition of two C-H bonds is inferred from alkene hydrogenation.", "Double-bond order decrease is not mislabeled as a new C-C sigma bond.", "No complete atom mapping, hydrogen source mapping, or coproduct balance."])
r["evidence_ids"] = ev
r["review"]["unresolved"] = ["Catalyst loading, H2 pressure, temperature, duration, solvent volume, workup, purification, yield basis, and selectivity method are unreported.", "R/S descriptors of the two centers are not printed on supplied pages."]
reactions.append(r)

# Scheme 3c, first arrow: Prins-type cyclization to 9.
r = set_common_reaction(st.new_reaction(f"{DOC}:scheme3c_prins_cyclization", DOC, ev), ev, ["derivatization"], {"page": 4, "scheme": "3c", "arrow": 1}, [], "Prins-type cyclization of 4aa with 3-nitrobenzaldehyde")
r["transformation_type"] = inferred("Prins-type cyclization", ev, "Normalized from the author's description.")
r["participants"] = [charge_part("p_4aa", "4aa", "reactant", ev), charge_part("p_8", "8", "reactant", ev), charge_part("p_bf3", "BF3OEt2", "reagent", ev), charge_part("p_thf", "THF", "solvent", ev), product_part("p_9", "9", ev)]
r["solvent_conditions"]["composition"] = f("THF", [molecules["THF"]["molecule_id"]], ev)
r["steps"] = [reaction_step(ev, [p["participant_id"] for p in r["participants"]], operation="Prins-type cyclization")]
r["outcomes"] = [outcome("yield", "yield", ["9"], ev, {"number": 75, "unit": "%"}, "bicyclic pyridinium salt 9", None, None, {"scheme": "3c first arrow"}), outcome("dr", "selectivity", ["9"], ev, "complete diastereoselectivity", "product 9", None, None, {"metric": "diastereoselectivity"})]
r["classification"] = classification(ev, [("C", "C", "formed"), ("C", "O", "formed")], limitations=["C-C and C-O ring closures are inferred from the Prins-type transformation drawing.", "Exact atom correspondences, proton transfers, and BF4 salt formation are unresolved."])
r["evidence_ids"] = ev
r["review"]["unresolved"] = ["Reagent amounts, temperature, duration, workup, purification, yield basis, and selectivity method are unreported.", "Complex product 9 connectivity and absolute descriptors require SI/expert structure review."]
reactions.append(r)

# Scheme 3c, second arrow: hydrolysis to 10.
r = set_common_reaction(st.new_reaction(f"{DOC}:scheme3c_hydrolysis", DOC, ev), ev, ["derivatization"], {"page": 4, "scheme": "3c", "arrow": 2}, [], "hydrolysis of bicyclic pyridinium salt 9")
r["transformation_type"] = inferred("hydrolysis", ev, "Normalized from the author's description.")
r["participants"] = [charge_part("p_9", "9", "reactant", ev), charge_part("p_k2co3", "K2CO3", "reagent", ev), product_part("p_10", "10", ev)]
r["steps"] = [reaction_step(ev, [p["participant_id"] for p in r["participants"]], operation="hydrolysis with aqueous K2CO3")]
r["outcomes"] = [outcome("yield", "yield", ["10"], ev, {"number": 85, "unit": "%"}, "product 10", None, None, {"scheme": "3c second arrow"}), outcome("er", "selectivity", ["10"], ev, "96:4 er", "product 10", None, None, {"metric": "er", "order": "major:minor"}), outcome("dr", "selectivity", ["10"], ev, "only diastereomer", "product 10", None, None, {"metric": "diastereoselectivity"})]
r["classification"] = classification(ev, coverage="unresolved", limitations=["The complex salt-to-pyridone atom mapping is not safely reconstructed from the compact main-text scheme.", "No bond-formation label is asserted until product 9/10 structures are expert-reviewed."])
r["evidence_ids"] = ev
r["review"]["unresolved"] = ["Aqueous K2CO3 amount/concentration, temperature, duration, workup, purification, yield basis, and selectivity methods are unreported.", "Structures 9 and 10 are source-identified but not encoded; bond classification remains unresolved."]
reactions.append(r)

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:page1_figure1", "kind": "background_examples", "evidence_ids": E[1], "note": "Figure 1 contains bioactive examples, not current-paper synthesis experiments."},
    {"context_id": f"{DOC}:page1_scheme1_prior", "kind": "prior_work", "evidence_ids": E[1], "note": "Scheme 1a and 1b summarize cited prior work; they are not counted as current-paper experiments."},
    {"context_id": f"{DOC}:page2_3_mechanistic_discussion", "kind": "author_interpretation", "evidence_ids": E[2] + E[3], "note": "Hydroxyl directing and hydrogen-bond explanations are author proposals, not separate experimental records."},
    {"context_id": f"{DOC}:si_notice", "kind": "supplementary_information_notice", "evidence_ids": E[4], "note": "The paper states that SI exists, but no SI file was explicitly supplied or associated in this workspace."},
    {"context_id": f"{DOC}:references_pages5_6", "kind": "bibliography", "evidence_ids": E[5] + E[6], "note": "Pages 5-6 are references and contain no current-paper experiments."},
]

st.write(RUN / "submitted.json", bundle)

diagnostics = {}
for key, m in molecules.items():
    smi = st.value(m["smiles"])
    if smi:
        diagnostics[key] = st.stereo_inspect(smi)
st.write(RUN / "raw_model/stereo_diagnostics.json", diagnostics)

inventory = {
    "document_id": DOC,
    "pages_inspected": [1, 2, 3, 4, 5, 6],
    "raw_counts": {
        "table1_rows": 18,
        "figure2_results": 12,
        "figure2_new_records_after_4aa_deduplication": 11,
        "figure3_results": 12,
        "scheme2_concrete_experiments": 8,
        "scheme3_concrete_experiments": 4,
        "specific_experiments": sum(x["record_type"] == "specific_experiment" for x in reactions),
        "generic_schemes": sum(x["record_type"] == "generic_scheme" for x in reactions),
        "proposed_mechanisms": sum(x["record_type"] == "proposed_mechanism" for x in reactions),
        "molecules": len(molecules),
    },
    "deduplication": {"Figure 2 product 4aa": f"{DOC}:table1_entry_09"},
    "source_limitations": ["No SI or CIF was supplied.", "Scheme 3a prints 'As Table 2' but no Table 2 is present in the supplied main-text pages.", "Many compact scope structures retain source labels and stereo descriptions without forced complete SMILES."],
}
st.write(RUN / "raw_model/extraction_inventory.json", inventory)

mentions = st.source_mentions(json.loads((RUN / "dataset.json").read_text()),
                              ["1a", "2a", "3a", "4aa", "1b", "6a", "9", "10"],
                              context=180)
st.write(RUN / "raw_model/source_mentions.json", mentions)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1,
            "state": "attempted",
            "regions": [
                {"region": "title, abstract, and graphical abstract", "finding": "Current-work generic aminoarylation family; graphical summary is not a separate experiment.", "evidence_ids": E[1]},
                {"region": "Figure 1", "finding": "Bioactive chiral N-substituted 2-pyridones are background examples, not current-paper experiments.", "evidence_ids": E[1]},
                {"region": "Scheme 1a and 1b", "finding": "Cited prior-work [3,3]-rearrangement and allylic substitution; retained as source_context, not current experiments.", "evidence_ids": E[1]},
                {"region": "Scheme 1c", "finding": "This-work generic transformation plus author-proposed allylpalladium/hydrogen-bonded mechanism; recorded separately as generic_scheme and proposed_mechanism.", "evidence_ids": E[1]},
            ],
            "record_ids": [rid_generic, f"{DOC}:proposed_mechanism"],
            "unresolved": ["Generic R/Ar/n structures are not complete atom-mapped molecules.", "The proposed intermediates were not isolated and are not treated as experiments."],
        },
        {
            "page": 2,
            "state": "attempted",
            "regions": [
                {"region": "Table 1 reaction header and ligand structures L1-L8", "finding": "Common optimization substrates/products and catalyst system; all explicit footnote quantities were structured.", "evidence_ids": E[2]},
                {"region": "Table 1 entries 1-18", "finding": "Eighteen independent optimization experiments; entries 4, 7, and 8 report NR and have no invented product.", "evidence_ids": E[2]},
                {"region": "Table 1 footnotes a-d", "finding": "Isolated yield, chiral-HPLC er, crude-1H-NMR b/l methods and standard charges/temperature/time were mapped to each applicable row.", "evidence_ids": E[2]},
                {"region": "scope discussion in lower text", "finding": "Narrative cross-reference to Figures 2 and 3; no extra experiment was created from prose repetitions.", "evidence_ids": E[2]},
            ],
            "record_ids": [f"{DOC}:table1_entry_{i:02d}" for i in range(1, 19)],
            "unresolved": ["NR analytical method and material balance are not reported.", "SI-only workup, purification, and detailed characterization are unavailable."],
        },
        {
            "page": 3,
            "state": "attempted",
            "regions": [
                {"region": "Figure 2 header, products 4aa-4la, and caption", "finding": "Twelve scope outcomes; 4aa exactly repeats Table 1 entry 9 and was merged as the same experiment, while 4ba-4la are eleven new records.", "evidence_ids": E[3]},
                {"region": "Figure 3 header, products 4ab-4kf, and caption", "finding": "Twelve independent substrate-scope experiments; condition A/B, isolated yields, er, and b/l methods were retained.", "evidence_ids": E[3]},
                {"region": "Scheme 2 general arrow", "finding": "Generic allene-variant aminoarylation scheme, separate from concrete rows.", "evidence_ids": E[3]},
                {"region": "Scheme 2 rows 1b-1e and products 6a-6e", "finding": "Eight concrete experiments: five successes and three NR controls (1c, 1d, 1e).", "evidence_ids": E[3]},
                {"region": "mechanistic prose below Figures/Scheme", "finding": "Hydroxyl directing and hydrogen-bond explanations are author interpretation; not separate experiments.", "evidence_ids": E[3]},
            ],
            "record_ids": [rid_generic, rid_scheme2, f"{DOC}:table1_entry_09"] +
                          [f"{DOC}:figure2_{x}" for x in fig2_products] +
                          [f"{DOC}:figure3_{x}" for x in fig3_products] +
                          [f"{DOC}:scheme2_{x}" for x in scheme2_success] +
                          [f"{DOC}:scheme2_{x}_NR" for x in ["1c", "1d", "1e"]],
            "unresolved": ["Scheme 2 does not state measurement methods for yield, er, or b/l.", "Many compact scope structures and absolute descriptors require expert review/SI; unsafe complete SMILES were not forced."],
        },
        {
            "page": 4,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 3a", "finding": "One 7 mmol scale-up experiment with substrate charges, 71%/1.2 g product, 96:4 er, and 79:21 b/l.", "evidence_ids": E[4]},
                {"region": "Scheme 3b", "finding": "One selective Pd/C-H2 hydrogenation experiment to 7; 95% and only-diastereomer result.", "evidence_ids": E[4]},
                {"region": "Scheme 3c first arrow", "finding": "One BF3 etherate/THF Prins-type cyclization to pyridinium salt 9; 75% with complete diastereoselectivity.", "evidence_ids": E[4]},
                {"region": "Scheme 3c second arrow", "finding": "One aqueous-K2CO3 hydrolysis to 10; 85%, 96:4 er, only diastereomer.", "evidence_ids": E[4]},
                {"region": "Associated Content and Accession Codes", "finding": "SI and CCDC records are announced, but no SI/CIF was explicitly supplied; no association was inferred.", "evidence_ids": E[4]},
            ],
            "record_ids": [f"{DOC}:scheme3a_scale_up", f"{DOC}:scheme3b_hydrogenation", f"{DOC}:scheme3c_prins_cyclization", f"{DOC}:scheme3c_hydrolysis"],
            "unresolved": ["Scheme 3a says 'As Table 2' while prose says 'standard conditions'; the mismatch is retained as a conflict.", "Scheme 3 omits most quantities, times, temperatures, workup/purification, and measurement methods.", "Structures 9 and 10 and the hydrolysis bond changes are not safely resolved from the compact main-text scheme."],
        },
        {
            "page": 5,
            "state": "attempted",
            "regions": [
                {"region": "Notes, Acknowledgments, and References 1-13", "finding": "No current-paper reaction experiment, table, scheme, mass spectrum, or characterization dataset; bibliography only after front matter.", "evidence_ids": E[5]},
            ],
            "record_ids": [],
            "unresolved": [],
        },
        {
            "page": 6,
            "state": "attempted",
            "regions": [
                {"region": "References 13-20", "finding": "Bibliography continuation only; cited reactions are not current-paper experiments.", "evidence_ids": E[6]},
            ],
            "record_ids": [],
            "unresolved": [],
        },
    ],
}
st.write(RUN / "coverage/pages_0001_0006.json", coverage)

tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    if task.get("document_id") == DOC:
        task["state"] = "attempted"
        task["reason"] = "Original page image and extracted text inspected; detailed region mapping is in coverage/pages_0001_0006.json."
st.write(RUN / "tasks.json", tasks)
print(json.dumps(inventory["raw_counts"], ensure_ascii=False, indent=2))
