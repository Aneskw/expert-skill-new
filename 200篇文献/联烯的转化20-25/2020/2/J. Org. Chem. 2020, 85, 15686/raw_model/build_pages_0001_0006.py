#!/usr/bin/env python3
"""Build the page 1--6 v3.1 extraction for doc_496d854c8fcb9457.

This is a source-derived machine submission.  It intentionally keeps every
record needs_review/partial and does not claim that RDKit parsing verifies a
drawn structure.
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
DOC = "doc_496d854c8fcb9457"
E = {p: {k: f"{DOC}:page{p}:{k}" for k in ("text", "image")} for p in range(1, 7)}
bundle = st.read(RUN / "raw_model/pages_0001_0006_template.json")
bundle["documents"][0]["title"] = "Metal-Free Intramolecular Aminophosphination of Allenes"


def q(number, unit, raw_text=None, **extra):
    value = {"number": number, "unit": unit}
    if raw_text is not None:
        value["raw_text"] = raw_text
    value.update(extra)
    return value


def f(raw=None, normalized=None, ev=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation,
                   status=status)


def inferred(normalized, ev, derivation, raw=None):
    return f(raw, normalized, ev, method="structure_derivation", derivation=derivation)


def set_source_fact(obj, raw, normalized, ev):
    obj.clear()
    obj.update(f(raw, normalized, ev))


def molecule(mid, label, ev, name=None, smiles=None, *, stereo=None, note=None,
             aliases=None):
    ev = list(dict.fromkeys(ev))
    m = st.new_molecule(f"{DOC}:{mid}", DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name, name, ev) if name else f()
    m["names"]["english"] = f(name, name, ev) if name else f()
    m["names"]["abbreviations"] = f(label, [label], ev)
    m["names"]["aliases"] = f(aliases, aliases, ev) if aliases else f()
    m["cas_numbers"] = f()
    m["identity_note"] = f(note, note, ev) if note else f()
    if smiles:
        deriv = "Candidate connectivity transcribed from the cited source drawing/name; RDKit parsing is not source verification."
        m["representation_type"] = "smiles"
        m["smiles"] = inferred(smiles, ev, deriv, raw={"source": "drawing/name", "candidate": smiles})
        diag = st.stereo_inspect(smiles)
        m["stereochemistry"]["identity_key"] = inferred(
            diag["stereo_graph_key"], ev, "RDKit stereo_graph_key calculated from the machine candidate SMILES.")
        m["stereochemistry"]["isomeric_smiles"] = inferred(
            diag["canonical_isomeric_smiles"], ev,
            "RDKit canonical isomeric SMILES calculated from the machine candidate; untagged source stereo remains unassigned.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f()
    stereo = stereo or {}
    # Every molecule receives explicit stereo states. Source-specified details
    # override structural diagnostics; absence of a source assignment is never
    # converted into racemic.
    ch = stereo.get("chirality_status", "unresolved" if not smiles else "unspecified")
    ge = stereo.get("geometry_status", "unresolved" if not smiles else "not_applicable")
    rp = stereo.get("representation_status", "unresolved" if not smiles else
                    ("partial" if ch in {"unspecified", "specified", "partially_specified"}
                     or ge in {"unspecified", "specified", "partially_specified"} else "not_applicable"))
    for key, value in (("chirality_status", ch), ("geometry_status", ge),
                       ("representation_status", rp)):
        m["stereochemistry"][key] = inferred(
            value, ev, "Status assigned conservatively from the source depiction/name and the candidate graph; missing tags are not treated as racemic.")
    for key in ("atom_centers", "double_bonds", "axial", "mixture", "assignment"):
        if key in stereo:
            m["stereochemistry"][key] = f(stereo[key], stereo[key], ev)
        else:
            m["stereochemistry"][key] = f()
    bundle["molecules"].append(m)
    return m["molecule_id"]


def template_molecule(mid, label, ev, fixed, variables, reason):
    ev = list(dict.fromkeys(ev))
    m = st.new_molecule(f"{DOC}:{mid}", DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(label, label, ev)
    m["names"]["english"] = f(label, label, ev)
    m["names"]["abbreviations"] = f(label, [label], ev)
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    m["representation_type"] = "generic_template"
    m["smiles"] = f()
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, ev),
        "variables": [{
            "label": v,
            "attachment_points": f("source drawing", "source drawing", ev),
            "definition": f(desc, desc, ev),
            "allowed_values": f(),
            "constraints": f(),
        } for v, desc in variables.items()],
        "instance_mapping": [],
        "reason": f(reason, reason, ev),
    }
    for key, value in (("chirality_status", "unresolved"),
                       ("geometry_status", "unresolved"),
                       ("representation_status", "not_encoded")):
        m["stereochemistry"][key] = inferred(
            value, ev, "Generic R-group drawing does not define one concrete stereoisomer.")
    m["identity_note"] = f(reason, reason, ev)
    bundle["molecules"].append(m)
    return m["molecule_id"]


M = {}


def add(key, *args, **kwargs):
    M[key] = molecule(key, *args, **kwargs)


def add_template(key, *args, **kwargs):
    M[key] = template_molecule(key, *args, **kwargs)


# Generic source drawings.
add_template("tmpl_1", "1 (generic sulfonamidoallene)", [E[1]["image"], E[2]["image"]],
             "Sulfonamide-tethered allene scaffold shown in Schemes 1/3",
             {"R1": "geminal tether substituent(s)", "R2/PG": "nitrogen protecting group"},
             "The figure contains R/PG variables and is not a concrete molecule.")
add_template("tmpl_3", "3 (generic pyrrolidine product)", [E[1]["image"], E[2]["image"]],
             "Alkenylphosphine-substituted pyrrolidine scaffold shown in Schemes 1/3",
             {"R1": "ring substituent(s)", "R2/PG": "nitrogen protecting group"},
             "The scope overview is a variable template, not one stereoisomer.")
add_template("tmpl_4", "4 (generic 2-allenylphenylmethanamine)", [E[2]["image"]],
             "2-Allenylphenylmethanamine scaffold shown in Scheme 4",
             {"R1": "alpha substituent", "R2": "aryl-side allene substituent", "R3/R4": "terminal allene substituents"},
             "Scheme 4 defines multiple variables and cannot be represented as one concrete SMILES.")
add_template("tmpl_5", "5 (generic isoindoline product)", [E[2]["image"]],
             "Alkenylphosphine-substituted fused N-heterocycle scaffold shown in Scheme 4",
             {"R1": "ring substituent", "R2/R3/R4": "alkene/ring substituents"},
             "Scheme 4 product is a variable structure and may include stereoisomers.")
add_template("tmpl_2", "2 (generic diarylphosphine oxide)", [E[3]["image"]],
             "H-P(=O)Ar2 shown in Scheme 5", {"Ar": "aryl group"},
             "Ar is varied in Scheme 5; the template is not a single substance.")
add_template("tmpl_6", "6 (generic phosphine-scope product)", [E[3]["image"]],
             "Isoindoline bearing C(=CH2)-PAr2 shown in Scheme 5", {"Ar": "aryl group"},
             "Ar is varied in Scheme 5; the template is not a single substance.")

# Concrete scope substrates and products.  Stereo states are deliberately
# conservative where the paper does not assign a configuration.
add("1a", "1a", [E[2]["image"], E[2]["text"]],
    "N-(hexa-4,5-dien-1-yl)-4-methylbenzenesulfonamide",
    "C=C=CCCCNS(=O)(=O)c1ccc(C)cc1",
    stereo={"chirality_status": "achiral", "geometry_status": "not_applicable",
            "representation_status": "not_applicable"})
add("1b", "1b", [E[2]["image"]], smiles="C=C=CCCCNS(=O)(=O)c1ccccc1",
    note="Scheme 3 structure candidate; full systematic substrate name is not printed in pages 1-6.",
    stereo={"chirality_status": "achiral", "geometry_status": "not_applicable", "representation_status": "not_applicable"})
add("1c", "1c", [E[2]["image"]], smiles="C=C=CCCCNS(=O)(=O)c1ccc([N+](=O)[O-])cc1",
    note="Nosyl-protected substrate structure candidate from Scheme 3.",
    stereo={"chirality_status": "achiral", "geometry_status": "not_applicable", "representation_status": "not_applicable"})
add("1d", "1d", [E[2]["image"]], smiles="C=C=CCC(C)(C)CNS(=O)(=O)c1ccc(C)cc1",
    note="Gem-dimethyl substrate structure candidate from Scheme 3.",
    stereo={"chirality_status": "achiral", "geometry_status": "not_applicable", "representation_status": "not_applicable"})
add("1e", "1e", [E[2]["image"]], smiles="C=C=CCC1(CNS(=O)(=O)c2ccc(C)cc2)CCCCC1",
    note="Spirocyclohexyl precursor structure candidate from Scheme 3.",
    stereo={"chirality_status": "achiral", "geometry_status": "not_applicable", "representation_status": "not_applicable"})
add("2a", "2a", [E[1]["image"], E[2]["image"], E[3]["text"]], "diphenylphosphine oxide",
    "O=P(c1ccccc1)c1ccccc1", stereo={"chirality_status": "achiral", "geometry_status": "not_applicable", "representation_status": "not_applicable"})

prod3 = {
    "3a": ("2-(1-(Diphenylphosphanyl)vinyl)-1-tosylpyrrolidine", "C=C(P(c1ccccc1)c1ccccc1)C1CCCN1S(=O)(=O)c1ccc(C)cc1"),
    "3ap": ("Diphenyl(1-(1-tosylpyrrolidin-2-yl)vinyl)phosphine oxide", "C=C(P(=O)(c1ccccc1)c1ccccc1)C1CCCN1S(=O)(=O)c1ccc(C)cc1"),
    "3b": ("2-(1-(Diphenylphosphanyl)vinyl)-1-(phenylsulfonyl)pyrrolidine", "C=C(P(c1ccccc1)c1ccccc1)C1CCCN1S(=O)(=O)c1ccccc1"),
    "3c": ("2-(1-(Diphenylphosphanyl)vinyl)-1-((4-nitrophenyl)sulfonyl)pyrrolidine", "C=C(P(c1ccccc1)c1ccccc1)C1CCCN1S(=O)(=O)c1ccc([N+](=O)[O-])cc1"),
    "3d": ("2-(1-(Diphenylphosphanyl)vinyl)-4,4-dimethyl-1-tosylpyrrolidine", "C=C(P(c1ccccc1)c1ccccc1)C1CC(C)(C)CN1S(=O)(=O)c1ccc(C)cc1"),
    "3e": ("3-(1-(Diphenylphosphanyl)vinyl)-2-tosyl-2-azaspiro[4.5]decane", "C=C(P(c1ccccc1)c1ccccc1)C1CC2(CCCCC2)CN1S(=O)(=O)c1ccc(C)cc1"),
}
for key, (name, smi) in prod3.items():
    lab = "3a′" if key == "3ap" else key
    page = 4
    add(key, lab, [E[2]["image"], E[page]["text"]], name, smi,
        stereo={"chirality_status": "unspecified", "geometry_status": "not_applicable", "representation_status": "partial"},
        note="The cyclized carbon is potentially stereogenic; pages 1-6 do not assign its absolute configuration.")

substrate_smiles = {
    "4a": "C=C=Cc1ccccc1CNS(=O)(=O)c1ccc(C)cc1",
    "4b": "C=C=Cc1cc(C)ccc1CNS(=O)(=O)c1ccc(C)cc1",
    "4c": "C=C=Cc1cc(OC)ccc1CNS(=O)(=O)c1ccc(C)cc1",
    "4d": "C=C=Cc1ccc(F)cc1CNS(=O)(=O)c1ccc(C)cc1",
    "4e": "C=C=Cc1ccc([N+](=O)[O-])cc1CNS(=O)(=O)c1ccc(C)cc1",
    "4f": "C=C=Cc1ccccc1C(C)NS(=O)(=O)c1ccc(C)cc1",
    "4g": "C=C=Cc1ccccc1C(CC)NS(=O)(=O)c1ccc(C)cc1",
    "4i": "CC=C=Cc1ccccc1CNS(=O)(=O)c1ccc(C)cc1",
    "4j": "CCCCCC=C=Cc1ccccc1CNS(=O)(=O)c1ccc(C)cc1",
    "4k": "CC(C)=C=Cc1ccccc1CNS(=O)(=O)c1ccc(C)cc1",
}
for key in [f"4{x}" for x in "abcdefg"]:
    stinfo = {"chirality_status": "unspecified" if key in {"4f", "4g"} else "achiral",
              "geometry_status": "not_applicable", "representation_status": "partial" if key in {"4f", "4g"} else "not_applicable"}
    add(key, key, [E[2]["image"]], smiles=substrate_smiles[key], stereo=stinfo,
        note="Scheme 4 substrate candidate; pages 1-6 do not print a full systematic name.")
add("4h", "4h", [E[2]["image"], E[3]["text"]], note="1,1-disubstituted allene; exact connectivity transcription remains unresolved from the compact source drawing.")
for key in ("4i", "4j"):
    add(key, key, [E[2]["image"], E[3]["text"]], smiles=substrate_smiles[key],
        stereo={"chirality_status": "unspecified", "geometry_status": "not_applicable", "representation_status": "not_encoded",
                "axial": [{"axis": "allene", "configuration": "not assigned in source"}]},
        note="Both allene termini have different substituents; axial configuration is not assigned in pages 1-6.")
add("4k", "4k", [E[2]["image"], E[3]["text"]], smiles=substrate_smiles["4k"],
    stereo={"chirality_status": "achiral", "geometry_status": "not_applicable", "representation_status": "not_applicable"},
    note="Trisubstituted allene candidate from Scheme 4.")
add("4l", "4l", [E[2]["image"], E[3]["text"]], "N-(2-(buta-2,3-dien-2-yl)benzyl)-4-methylbenzenesulfonamide",
    note="The compact drawing/regiochemical relation to 5l was not encoded as a concrete SMILES without SI.")

prod5 = {
    "5a": ("1-(1-(Diphenylphosphanyl)vinyl)-2-tosylisoindoline", "C=C(P(c1ccccc1)c1ccccc1)C1N(S(=O)(=O)c2ccc(C)cc2)Cc2ccccc21"),
    "5b": ("1-(1-(Diphenylphosphanyl)vinyl)-5-methyl-2-tosylisoindoline", "C=C(P(c1ccccc1)c1ccccc1)C1N(S(=O)(=O)c2ccc(C)cc2)Cc2ccc(C)cc21"),
    "5c": ("1-(1-(Diphenylphosphanyl)vinyl)-5-methoxy-2-tosylisoindoline", "C=C(P(c1ccccc1)c1ccccc1)C1N(S(=O)(=O)c2ccc(C)cc2)Cc2ccc(OC)cc21"),
    "5d": ("1-(1-(Diphenylphosphanyl)vinyl)-6-fluoro-2-tosylisoindoline", "C=C(P(c1ccccc1)c1ccccc1)C1N(S(=O)(=O)c2ccc(C)cc2)Cc2cc(F)ccc21"),
    "5e": ("1-(1-(Diphenylphosphanyl)vinyl)-6-nitro-2-tosylisoindoline", "C=C(P(c1ccccc1)c1ccccc1)C1N(S(=O)(=O)c2ccc(C)cc2)Cc2cc([N+](=O)[O-])ccc21"),
    "5f": ("(1R,3S)-1-(1-(Diphenylphosphanyl)vinyl)-3-methyl-2-tosylisoindoline", "C=C(P(c1ccccc1)c1ccccc1)C1N(S(=O)(=O)c2ccc(C)cc2)C(C)c2ccccc21"),
    "5g": ("(1R,3S)-1-(1-(Diphenylphosphanyl)vinyl)-3-ethyl-2-tosylisoindoline", "C=C(P(c1ccccc1)c1ccccc1)C1N(S(=O)(=O)c2ccc(C)cc2)C(CC)c2ccccc21"),
    "5h": ("1-(1-(Diphenylphosphanyl)vinyl)-1-methyl-2-tosylisoindoline", "C=C(P(c1ccccc1)c1ccccc1)C1(C)N(S(=O)(=O)c2ccc(C)cc2)Cc2ccccc21"),
    "5i": ("(E)-1-(1-(Diphenylphosphanyl)prop-1-en-1-yl)-2-tosylisoindoline", "C/C=C(/C1N(S(=O)(=O)c2ccc(C)cc2)Cc2ccccc21)P(c1ccccc1)c1ccccc1"),
    "5j": ("(E)-1-(1-(Diphenylphosphanyl)hept-1-en-1-yl)-2-tosylisoindoline", "CCCCC/C=C(/C1N(S(=O)(=O)c2ccc(C)cc2)Cc2ccccc21)P(c1ccccc1)c1ccccc1"),
    "5k": ("1-(2-Methylprop-1-en-1-yl)-2-tosylisoindoline", "CC(C)=CC1N(S(=O)(=O)c2ccc(C)cc2)Cc2ccccc21"),
    "5l": ("1-(1-(Diphenylphosphanyl)vinyl)-2-tosyl-1,2,3,4-tetrahydroisoquinoline", "C=C(P(c1ccccc1)c1ccccc1)C1N(S(=O)(=O)c2ccc(C)cc2)CCc2ccccc21"),
}
for idx, key in enumerate([f"5{x}" for x in "abcdefghijkl"]):
    page = 4 if key in {"5a", "5b"} else 5 if key not in {"5l"} else 6
    name, smi = prod5[key]
    stere = {"chirality_status": "unspecified", "geometry_status": "not_applicable", "representation_status": "partial"}
    note = "The ring carbon may be stereogenic; no absolute configuration is assigned in pages 1-6."
    if key in {"5f", "5g"}:
        stere.update(chirality_status="specified", representation_status="not_encoded",
                      assignment={"source_name": name, "relative_configuration": "anti" if key == "5g" else "source name gives (1R,3S)"})
        note = "The source name gives (1R,3S); the machine connectivity SMILES does not encode the two centers. anti-5g is assigned by single-crystal X-ray diffraction."
    if key in {"5i", "5j"}:
        stere.update(geometry_status="specified", representation_status="partial",
                      double_bonds=[{"alkene": "P-substituted exocyclic alkene", "configuration": "E"}])
        note = "E geometry is source-labelled; 5i is stated to be established by NOESY. The ring stereocenter remains unspecified."
    if key == "5k":
        note = "Hydroamination product; no C-P bond is present. The ring stereocenter is not assigned."
    add(key, key, [E[2]["image"], E[page]["text"]], name, smi, stereo=stere, note=note)

add("2b", "2 (di-4-tolyl)", [E[3]["image"], E[3]["text"]], "di-4-tolylphosphine oxide",
    "O=P(c1ccc(C)cc1)c1ccc(C)cc1", stereo={"chirality_status": "achiral", "geometry_status": "not_applicable", "representation_status": "not_applicable"})
add("2c", "2 (di([1,1′-biphenyl]-4-yl))", [E[3]["image"], E[3]["text"]], "di([1,1′-biphenyl]-4-yl)phosphine oxide",
    "O=P(c1ccc(-c2ccccc2)cc1)c1ccc(-c2ccccc2)cc1", stereo={"chirality_status": "achiral", "geometry_status": "not_applicable", "representation_status": "not_applicable"})
add("2d", "2 (di(naphthalen-1-yl))", [E[3]["image"]], "di(naphthalen-1-yl)phosphine oxide",
    "O=P(c1cccc2ccccc12)c1cccc2ccccc12", stereo={"chirality_status": "unspecified", "geometry_status": "not_applicable", "representation_status": "partial"},
    note="Possible phosphorus configurational issues are not assigned; source specifies only the aryl connectivity.")
add("2e", "2 (di(naphthalen-2-yl))", [E[3]["image"]], "di(naphthalen-2-yl)phosphine oxide",
    "O=P(c1ccc2ccccc2c1)c1ccc2ccccc2c1", stereo={"chirality_status": "unspecified", "geometry_status": "not_applicable", "representation_status": "partial"})

prod6 = {
    "6a": ("1-(1-(Di-p-tolylphosphanyl)vinyl)-2-tosylisoindoline", "C=C(P(c1ccc(C)cc1)c1ccc(C)cc1)C1N(S(=O)(=O)c2ccc(C)cc2)Cc2ccccc21"),
    "6b": ("1-(1-(Di([1,1′-biphenyl]-4-yl)phosphanyl)vinyl)-2-tosylisoindoline", "C=C(P(c1ccc(-c2ccccc2)cc1)c1ccc(-c2ccccc2)cc1)C1N(S(=O)(=O)c2ccc(C)cc2)Cc2ccccc21"),
    "6c": ("1-(1-(Di(naphthalen-1-yl)phosphanyl)vinyl)-2-tosylisoindoline", "C=C(P(c1cccc2ccccc12)c1cccc2ccccc12)C1N(S(=O)(=O)c2ccc(C)cc2)Cc2ccccc21"),
    "6d": ("1-(1-(Di(naphthalen-2-yl)phosphanyl)vinyl)-2-tosylisoindoline", "C=C(P(c1ccc2ccccc2c1)c1ccc2ccccc2c1)C1N(S(=O)(=O)c2ccc(C)cc2)Cc2ccccc21"),
}
for key, (name, smi) in prod6.items():
    add(key, key, [E[3]["image"], E[6]["text"]], name, smi,
        stereo={"chirality_status": "unspecified", "geometry_status": "not_applicable", "representation_status": "partial"},
        note="Potential ring stereocenter is not assigned in pages 1-6.")

# Trace controls and common process materials.
add("phosphonate_control", "HP(O)Ph(OEt)", [E[3]["text"]], "ethyl phenylphosphinate",
    "CCOP(=O)c1ccccc1", stereo={"chirality_status": "unspecified", "geometry_status": "not_applicable", "representation_status": "partial"},
    note="Name/formula is transcribed as HP(O)Ph(OEt); candidate SMILES retains an implicit P-H.")
add("oxaphosphinine_control", "dibenzo[c,e][1,2]oxaphosphinine 6-oxide", [E[3]["text"]],
    "dibenzo[c,e][1,2]oxaphosphinine 6-oxide", note="Compact heterocycle structure is not drawn in the supplied pages and is left unresolved.")
add("trace_product_1", "trace aminophosphination product (HP(O)Ph(OEt))", [E[3]["text"]], note="Trace product is mentioned but not structurally identified.")
add("trace_product_2", "trace aminophosphination product (oxaphosphinine source)", [E[3]["text"]], note="Trace product is mentioned but not structurally identified.")

common = {
    "tf2o": ("Tf2O", "triflic anhydride", "O=S(=O)(OS(=O)(=O)C(F)(F)F)C(F)(F)F"),
    "lutidine": ("2,6-lutidine", "2,6-lutidine", "Cc1cccc(C)n1"),
    "dce": ("DCE", "1,2-dichloroethane", "ClCCCl"),
    "n2": ("N2", "nitrogen", "N#N"),
    "h2o2": ("H2O2", "hydrogen peroxide", "OO"),
    "nahco3": ("sat. NaHCO3", "saturated aqueous sodium bicarbonate", "O=C([O-])O.[Na+]"),
    "dcm": ("CH2Cl2", "dichloromethane", "ClCCl"),
    "na2so4": ("anhydrous Na2SO4", "anhydrous sodium sulfate", "O=S(=O)([O-])[O-].[Na+].[Na+]"),
    "etoac": ("EtOAc", "ethyl acetate", "CCOC(C)=O"),
}
for key, (label, name, smi) in common.items():
    ev = [E[3]["text"], E[4]["text"]] if key not in {"h2o2"} else [E[2]["image"], E[2]["text"]]
    add(key, label, ev, name, smi,
        stereo={"chirality_status": "achiral", "geometry_status": "not_applicable", "representation_status": "not_applicable"})
add_template("petroleum_ether", "PE", [E[4]["text"], E[5]["text"], E[6]["text"]],
             "petroleum ether, composition not reported", {"composition": "unreported hydrocarbon mixture"},
             "A named solvent mixture cannot be assigned a unique molecular SMILES.")
add_template("silica", "silica gel", [E[3]["text"], E[4]["text"]],
             "silica gel (200-300 mesh)", {"surface": "heterogeneous material"},
             "Silica gel is a material rather than a single molecular species.")


def set_charge(p, key, raw, normalized, ev, *, reference=None, derived=False):
    if derived:
        p["charges"][key] = inferred(normalized, ev, "Calculated from the explicitly reported charge ratio.", raw=raw)
    else:
        p["charges"][key] = f(raw, normalized, ev)
    if reference:
        p["equivalence_reference"] = f(reference, reference, ev)


def participant(r, pid, mkey, role, ev, stage="reaction"):
    p = st.new_participant(pid, M[mkey], role, ev, stage=stage)
    r["participants"].append(p)
    return p


def change(change_id, a, b, kind="formed", before=0, after=1):
    return {"change_id": change_id, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
            "before_order": before, "after_order": after, "change_type": kind}


def set_classification(r, pairs, ev, note):
    changes = [change(f"b{i+1}", a, b) for i, (a, b) in enumerate(pairs)]
    labels = []
    for i, (a, b) in enumerate(pairs):
        lab = st.bond_label(a, b, 1)
        existing = next((x for x in labels if x["code"] == lab["code"]), None)
        if existing:
            existing["supporting_change_ids"].append(f"b{i+1}")
        else:
            labels.append({**lab, "supporting_change_ids": [f"b{i+1}"]})
    r["classification"]["mapped_reaction_smiles"] = f()
    r["classification"]["bond_changes"] = [inferred(c, ev, note) for c in changes]
    r["classification"]["labels"] = inferred(labels, ev, note)
    r["classification"]["coverage"] = inferred("partial", ev, note)
    r["classification"]["limitations"] = [
        "No complete atom mapping was asserted.",
        "Labels represent the visually supported net newly formed bonds only; omitted coproducts and proton transfers remain unresolved.",
    ]
    r["classification"]["model_votes"] = []


def base_reaction(rid, record_type, ev, location, reaction_name, transform, categories=None):
    ev = list(dict.fromkeys(ev))
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, ev) if categories is not None else f()
    r["location"] = f(location, location, ev)
    r["reaction_name"] = f(reaction_name, reaction_name, ev) if reaction_name else f()
    r["transformation_type"] = inferred(transform, ev, "Machine normalization of the source reaction description.") if transform else f()
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["Machine extraction requires expert review; schema validation is not chemical verification."]}
    r["completeness"] = "partial"
    return r


def add_overview(rid, rtype, page, location, title, raw_note, participants=(), pairs=()):
    ev = [E[page]["image"], E[page]["text"]]
    r = base_reaction(rid, rtype, ev, location, title, None)
    for i, (mkey, role) in enumerate(participants):
        participant(r, f"p{i+1}", mkey, role, ev)
    r["notes"] = [f(raw_note, raw_note, ev)]
    if pairs:
        set_classification(r, pairs, ev, "Net bond classes inferred from the generic arrow; variable identities and atom mapping remain unresolved.")
    else:
        r["classification"]["coverage"] = inferred("unresolved", ev, "Generic/prior-literature arrow was inventoried without asserting a complete bond classification.")
    r["review"]["unresolved"].append("Generic or prior-literature scheme is not a concrete experiment from this paper.")
    bundle["reactions"].append(r)
    return r["reaction_id"]


# Page-1 prior literature arrows and the paper's generic overview.
for suffix, desc in {
    "scheme1a_prior_cu": "Cu-catalyzed cyclization of Ts-protected o-alkenylaniline with diphenylphosphine oxide.",
    "scheme1b_prior_oxime": "Prior oxime-phosphonate alkene aminophosphonation example.",
    "scheme1c_prior_amino": "Prior Cu/Ag or Fe aminophosphinoylation of styrene/aniline systems.",
    "scheme1d_prior_azido": "Prior Mn-mediated phosphonation-azidation of alkenes.",
    "scheme1e_prior_alkyne": "Prior Tf2O/2,6-lutidine electrophilic phosphination/cyclization of alkynes.",
}.items():
    add_overview(suffix, "literature_example", 1, f"PDF page 1, Scheme 1 ({suffix.split('_')[0][-1:] or suffix})", "Aminophosphination/Aminophosphinoylation of C-C Multiple Bonds", desc)
add_overview("scheme1_this_work", "generic_scheme", 1, "PDF page 1, Scheme 1, this work", "Metal-free intramolecular aminophosphination of allenes",
             "Generic overview: sulfonamidoallene plus diarylphosphine oxide with Tf2O/2,6-lutidine gives an alkenylphosphine N-heterocycle.",
             [("tmpl_1", "reactant"), ("tmpl_2", "reactant"), ("tmpl_3", "product")], [("C", "N"), ("C", "P")])
add_overview("scheme3_generic", "generic_scheme", 2, "PDF page 2, Scheme 3 generic arrow", "Scope of Intramolecular Aminophosphination of Allenes",
             "Variable sulfonamidoallene 1 and 2a give pyrrolidine products 3 under Tf2O/2,6-lutidine in DCE/N2 at room temperature for 2 h.",
             [("tmpl_1", "reactant"), ("2a", "reactant"), ("tmpl_3", "product")], [("C", "N"), ("C", "P")])
add_overview("scheme4_generic", "generic_scheme", 2, "PDF page 2, Scheme 4 generic arrow", "Scope of Intramolecular Aminophosphination of 2-Allenylphenylmethanamine",
             "Variable 2-allenylphenylmethanamine 4 and 2a give fused N-heterocycles 5 under the general conditions.",
             [("tmpl_4", "reactant"), ("2a", "reactant"), ("tmpl_5", "product")], [("C", "N"), ("C", "P")])
add_overview("scheme5_generic", "generic_scheme", 3, "PDF page 3, Scheme 5 generic arrow", "Scope of Phosphine Oxides 2",
             "4a and generic diarylphosphine oxide 2 give product 6 under the general conditions.",
             [("4a", "reactant"), ("tmpl_2", "reactant"), ("tmpl_6", "product")], [("C", "N"), ("C", "P")])


GENERAL_TEXT = (
    "In a Schlenk tube, allene substrates (0.2 mmol, 1.0 equiv) and phosphine oxides 2 "
    "(0.4 mmol, 2.0 equiv) were placed with a magnetic stir bar under a N2 atmosphere. "
    "1,2-Dichloroethane (2 mL, 0.1 M) and 2,6-lutidine (47 μL, 0.4 mmol) were added by a syringe sequentially. "
    "After homogenization, Tf2O (67 μL, 0.4 mmol) was added by syringe and the mixture was stirred 2 h at room temperature; "
    "it was poured into saturated NaHCO3, extracted with CH2Cl2, dried over anhydrous Na2SO4, concentrated, and purified by short silica-gel chromatography."
)


def scope_reaction(rid, substrate, phosphine, product, product_page, mass_mg, yield_pct,
                   ratio="8:1", categories=("substrate_scope",), extra=None):
    scheme_page = 2 if rid.startswith(("3", "5")) else 3
    ev = [E[scheme_page]["image"], E[3]["text"], E[4]["text"], E[product_page]["text"]]
    # Scheme/text discussion on page 3 is relevant for 5g-l and phosphine controls.
    if rid in {"5g", "5h", "5i", "5j", "5k", "5l"} and E[3]["text"] not in ev:
        ev.append(E[3]["text"])
    ev = list(dict.fromkeys(ev))
    r = base_reaction(rid, "specific_experiment", ev,
                      f"Scheme {'3' if rid.startswith('3') else '4' if rid.startswith('5') else '5'}; product paragraph on PDF page {product_page}",
                      "intramolecular aminophosphination of allenes" if rid != "5k" else "intramolecular hydroamination of an allene",
                      "allene aminophosphination/cyclization" if rid != "5k" else "allene hydroamination/cyclization",
                      list(categories))
    ps = participant(r, "substrate", substrate, "reactant", ev)
    pp = participant(r, "phosphine", phosphine, "reactant" if rid != "5k" else "reagent", ev)
    pl = participant(r, "lutidine", "lutidine", "reagent", ev)
    pt = participant(r, "tf2o", "tf2o", "reagent", ev)
    pd = participant(r, "dce", "dce", "solvent", ev)
    pn = participant(r, "nitrogen", "n2", "additive", ev)
    po = participant(r, "product", product, "product", ev)
    pw1 = participant(r, "nahco3", "nahco3", "workup_agent", ev, "workup")
    pw2 = participant(r, "dcm", "dcm", "workup_agent", ev, "workup")
    pw3 = participant(r, "na2so4", "na2so4", "workup_agent", ev, "workup")
    ppe = participant(r, "pe", "petroleum_ether", "purification_agent", ev, "purification")
    pea = participant(r, "etoac", "etoac", "purification_agent", ev, "purification")
    psi = participant(r, "silica", "silica", "purification_agent", ev, "purification")
    set_charge(ps, "amount", "0.2 mmol", q(0.2, "mmol"), [E[3]["text"]])
    set_charge(ps, "equivalents", "1.0 equiv", q(1.0, "equiv"), [E[3]["text"]], reference=M[substrate])
    set_charge(pp, "amount", "0.4 mmol", q(0.4, "mmol"), [E[3]["text"]])
    set_charge(pp, "equivalents", "2.0 equiv", q(2.0, "equiv"), [E[3]["text"]], reference=M[substrate])
    set_charge(pl, "amount", "0.4 mmol", q(0.4, "mmol"), [E[3]["text"]])
    set_charge(pl, "volume", "47 μL", q(47, "uL"), [E[3]["text"]])
    set_charge(pt, "amount", "0.4 mmol", q(0.4, "mmol"), [E[4]["text"]])
    set_charge(pt, "volume", "67 μL", q(67, "uL"), [E[4]["text"]])
    set_charge(pd, "volume", "2 mL", q(2, "mL"), [E[3]["text"]])
    set_charge(po, "mass", f"{mass_mg} mg", q(mass_mg, "mg"), [E[product_page]["text"]])
    r["reaction_scale"] = f("allene substrate 0.2 mmol", {"reference_molecule_id": M[substrate], "amount": q(0.2, "mmol")}, [E[3]["text"]])
    r["solvent_conditions"] = {
        "composition": f("1,2-dichloroethane", [M["dce"]], [E[3]["text"]]),
        "mixture_ratio": f(),
        "total_volume": f("2 mL", q(2, "mL"), [E[3]["text"]]),
        "concentration": f("0.1 M", q(0.1, "M"), [E[3]["text"]]),
        "concentration_basis": f("allene substrate", M[substrate], [E[3]["text"]]),
    }
    s1 = st.new_step("charge_under_n2", ev)
    s1["participant_ids"] = [ps["participant_id"], pp["participant_id"], pn["participant_id"]]
    s1["operation"] = f("placed with a magnetic stir bar under a N2 atmosphere", "charge substrate and phosphine oxide under N2", [E[3]["text"]])
    s1["atmosphere"] = f("N2", "N2", [E[3]["text"]])
    s2 = st.new_step("add_solvent_and_base", ev)
    s2["participant_ids"] = [pd["participant_id"], pl["participant_id"]]
    s2["operation"] = f("added by a syringe sequentially; stirred until homogeneous", "sequential syringe addition then homogenization", [E[3]["text"], E[4]["text"]])
    s2["atmosphere"] = f("N2", "N2", [E[3]["text"]])
    s3 = st.new_step("add_tf2o_and_hold", ev)
    s3["participant_ids"] = [pt["participant_id"]]
    s3["operation"] = f("Tf2O was then added by a syringe; stirred", "syringe addition and hold", [E[4]["text"]])
    s3["temperature"] = f("room temperature", "room temperature", [E[4]["text"]])
    s3["time"] = f("2 h", q(2, "h"), [E[4]["text"]])
    s3["atmosphere"] = f("N2", "N2", [E[3]["text"]])
    s4 = st.new_step("workup", ev, stage="workup")
    s4["participant_ids"] = [pw1["participant_id"], pw2["participant_id"], pw3["participant_id"]]
    s4["operation"] = f("poured into sat. NaHCO3; extracted with CH2Cl2; combined organic layer dried over anhydrous Na2SO4; volatiles removed under reduced pressure",
                        "quench, extract, combine, dry, and concentrate", [E[4]["text"]])
    s5 = st.new_step("purification", ev, stage="purification")
    s5["participant_ids"] = [ppe["participant_id"], pea["participant_id"], psi["participant_id"]]
    s5["operation"] = f(f"short column chromatography (PE/EtOAc = {ratio})", {"method": "short silica-gel column chromatography", "eluent": "PE/EtOAc", "ratio": ratio}, [E[product_page]["text"]])
    r["steps"] = [s1, s2, s3, s4, s5]
    r["procedure"]["original_text"] = f(GENERAL_TEXT, GENERAL_TEXT, [E[3]["text"], E[4]["text"]])
    r["procedure"]["general_procedure_reference"] = f("prepared according to the general procedure", "General Procedures", [E[product_page]["text"]])
    r["procedure"]["addition_order"] = f("allene + phosphine oxide; DCE; 2,6-lutidine; homogenize; Tf2O", [M[substrate], M[phosphine], M["dce"], M["lutidine"], M["tf2o"]], [E[3]["text"], E[4]["text"]])
    oy = st.new_outcome("isolated_yield", "yield", [M[product]], ev)
    oy["product_role"] = f("desired product" if rid != "5k" else "hydroamination product", "main_product" if rid != "5k" else "side_product", [E[product_page]["text"]])
    oy["value"] = f(f"{yield_pct}% yield", q(yield_pct, "%"), [E[product_page]["text"]])
    oy["yield_type"] = f("isolated", "isolated", [E[product_page]["text"]])
    oy["measurement_basis"] = f("isolated after short column chromatography", "isolated after chromatography", [E[product_page]["text"]])
    om = st.new_outcome("product_amount", "product_amount", [M[product]], ev)
    om["product_role"] = f("isolated product", "main_product" if rid != "5k" else "side_product", [E[product_page]["text"]])
    om["value"] = f(f"{mass_mg} mg", q(mass_mg, "mg"), [E[product_page]["text"]])
    om["measurement_basis"] = f("isolated after short column chromatography", "isolated mass", [E[product_page]["text"]])
    r["outcomes"] = [oy, om]
    if extra:
        extra(r, ev)
    set_classification(r, [("C", "N"), ("C", "H")] if rid == "5k" else [("C", "N"), ("C", "P")],
                       [E[scheme_page]["image"]],
                       "Net newly formed bonds inferred from the source substrate/product drawings; exact atom mapping and omitted byproducts are unresolved.")
    r["review"]["unresolved"] += [
        "Candidate reactant/product SMILES require expert comparison to the source drawing and unavailable SI.",
        "Exact atom mapping, all proton/oxygen transfers, and coproduct identities are unresolved.",
        "Supporting Information was mentioned by the paper but was not supplied or associated in this dataset.",
    ]
    bundle["reactions"].append(r)
    return r


def h2o2_extra(r, ev):
    ph = participant(r, "h2o2", "h2o2", "reagent", [E[2]["image"]])
    r["steps"].insert(4, st.new_step("h2o2_treatment", [E[2]["image"]]))
    r["steps"][4]["participant_ids"] = [ph["participant_id"]]
    r["steps"][4]["operation"] = f("The reaction mixture was treated with H2O2.", "post-reaction H2O2 treatment", [E[2]["image"], E[2]["text"]])
    r["review"]["unresolved"].append("H2O2 amount, temperature, and treatment time are not reported in pages 1-6.")


def dr_extra(product):
    def apply(r, ev):
        o = st.new_outcome("diastereomeric_ratio", "selectivity", [M[product]], ev)
        o["value"] = f("dr >10:1", {"metric": "dr", "relation": ">", "ratio": [10, 1], "major_minor_assignment": "not fully specified"}, [E[2]["image"]])
        o["measurement_basis"] = f("calculated through the NMR spectrum of the crude product mixture", "crude-mixture NMR", [E[2]["image"]])
        o["scope"] = f("dr; major/minor configurations not both identified", {"metric": "dr", "order": "major:minor", "target": M[product]}, [E[2]["image"]])
        r["outcomes"].append(o)
    return apply


def ez_extra(product):
    def apply(r, ev):
        o = st.new_outcome("E_Z_ratio", "selectivity", [M[product]], ev)
        o["value"] = f("E:Z >10:1", {"metric": "E/Z", "relation": ">", "ratio": [10, 1], "order": ["E", "Z"]}, [E[2]["image"]])
        o["measurement_basis"] = f("calculated through the NMR spectrum of the crude product mixture", "crude-mixture NMR", [E[2]["image"]])
        o["scope"] = f("E:Z", {"metric": "E/Z", "order": ["E", "Z"], "target": M[product]}, [E[2]["image"]])
        r["outcomes"].append(o)
    return apply


scope_reaction("3a", "1a", "2a", "3a", 4, 75, 86)
scope_reaction("3ap", "1a", "2a", "3ap", 4, 65, 72, ratio="2:1", categories=("substrate_scope", "derivatization"), extra=h2o2_extra)
for rid, sub, mass, yp in [("3b", "1b", 60, 71), ("3c", "1c", 60, 65), ("3d", "1d", 70, 76), ("3e", "1e", 80, 79)]:
    scope_reaction(rid, sub, "2a", rid, 4, mass, yp)

five_data = [
    ("5a", "4a", 4, 90, 93, None), ("5b", "4b", 4, 90, 90, None),
    ("5c", "4c", 5, 95, 93, None), ("5d", "4d", 5, 85, 85, None),
    ("5e", "4e", 5, 75, 71, None), ("5f", "4f", 5, 85, 85, dr_extra("5f")),
    ("5g", "4g", 5, 90, 88, dr_extra("5g")), ("5h", "4h", 5, 80, 80, None),
    ("5i", "4i", 5, 85, 85, ez_extra("5i")), ("5j", "4j", 5, 90, 81, ez_extra("5j")),
    ("5k", "4k", 5, 57, 87, None), ("5l", "4l", 6, 90, 90, None),
]
for rid, sub, page, mass, yp, extra in five_data:
    scope_reaction(rid, sub, "2a", rid, page, mass, yp, ratio="10:1" if rid == "5k" else "8:1", extra=extra)

for rid, phosphine, mass, yp in [("6a", "2b", 70, 69), ("6b", "2c", 100, 79), ("6c", "2d", 108, 92), ("6d", "2e", 90, 77)]:
    scope_reaction(rid, "4a", phosphine, rid, 6, mass, yp)


# Gram-scale 5a is a distinct experiment, not the 0.2 mmol scope preparation.
ev = [E[4]["text"]]
r = base_reaction("5a_gram_scale", "specific_experiment", ev, "PDF page 4, Gram-Scale Preparation of 5a",
                  "Gram-Scale Preparation of 5a", "allene aminophosphination/cyclization", ["scale_up"])
parts = {}
for pid, key, role, stage in [
    ("substrate", "4a", "reactant", "reaction"), ("phosphine", "2a", "reactant", "reaction"),
    ("lutidine", "lutidine", "reagent", "reaction"), ("tf2o", "tf2o", "reagent", "reaction"),
    ("dce", "dce", "solvent", "reaction"), ("nitrogen", "n2", "additive", "reaction"),
    ("product", "5a", "product", "reaction"), ("nahco3", "nahco3", "workup_agent", "workup"),
    ("dcm", "dcm", "workup_agent", "workup"), ("na2so4", "na2so4", "workup_agent", "workup"),
    ("silica", "silica", "purification_agent", "purification")]:
    parts[pid] = participant(r, pid, key, role, ev, stage)
for pid, vals in {
    "substrate": [("mass", "1.5 g", q(1.5, "g")), ("amount", "5.0 mmol", q(5.0, "mmol"))],
    "phosphine": [("mass", "2.0 g", q(2.0, "g")), ("amount", "9.9 mmol", q(9.9, "mmol"))],
    "lutidine": [("volume", "1.2 mL", q(1.2, "mL")), ("amount", "10.3 mmol", q(10.3, "mmol"))],
    "tf2o": [("volume", "1.7 mL", q(1.7, "mL")), ("amount", "10.1 mmol", q(10.1, "mmol"))],
    "dce": [("volume", "50 mL", q(50, "mL"))],
    "product": [("mass", "2.1 g", q(2.1, "g"))],
}.items():
    for key, raw, norm in vals:
        set_charge(parts[pid], key, raw, norm, ev)
r["reaction_scale"] = f("4a (5.0 mmol)", {"reference_molecule_id": M["4a"], "amount": q(5.0, "mmol")}, ev)
r["solvent_conditions"] = {"composition": f("1,2-dichloroethane", [M["dce"]], ev), "mixture_ratio": f(),
                           "total_volume": f("50 mL", q(50, "mL"), ev), "concentration": f("0.1 M", q(0.1, "M"), ev),
                           "concentration_basis": f("4a", M["4a"], ev)}
for sid, stage, pids, op, temp, time in [
    ("charge_under_n2", "reaction", ["substrate", "phosphine", "nitrogen"], "placed with stir bar under N2", None, None),
    ("add_solvent_and_base", "reaction", ["dce", "lutidine"], "added sequentially by syringe and homogenized", None, None),
    ("add_tf2o_and_hold", "reaction", ["tf2o"], "Tf2O added by syringe; stirred", "room temperature", q(2, "h")),
    ("workup", "workup", ["nahco3", "dcm", "na2so4"], "poured into sat. NaHCO3, extracted with CH2Cl2, dried over anhydrous Na2SO4, concentrated", None, None),
    ("purification", "purification", ["silica"], "short column chromatography on silica gel", None, None),
]:
    s = st.new_step(sid, ev, stage=stage); s["participant_ids"] = [parts[x]["participant_id"] for x in pids]
    s["operation"] = f(op, op, ev)
    if temp: s["temperature"] = f(temp, temp, ev)
    if time: s["time"] = f("2 h", time, ev)
    if stage == "reaction": s["atmosphere"] = f("N2", "N2", ev)
    r["steps"].append(s)
r["procedure"]["original_text"] = f("Gram-Scale Preparation of 5a paragraph", "Full gram-scale procedure on PDF page 4", ev)
r["procedure"]["addition_order"] = f("4a + 2a; DCE; 2,6-lutidine; homogenize; Tf2O", [M["4a"], M["2a"], M["dce"], M["lutidine"], M["tf2o"]], ev)
oy = st.new_outcome("isolated_yield", "yield", [M["5a"]], ev); oy["value"] = f("88% yield", q(88, "%"), ev); oy["yield_type"] = f("isolated", "isolated", ev); oy["measurement_basis"] = f("white solid after short column chromatography", "isolated after chromatography", ev); oy["product_role"] = f("desired product 5a", "main_product", ev)
om = st.new_outcome("product_amount", "product_amount", [M["5a"]], ev); om["value"] = f("2.1 g", q(2.1, "g"), ev); om["measurement_basis"] = f("white solid after short column chromatography", "isolated mass", ev); om["product_role"] = f("desired product 5a", "main_product", ev)
r["outcomes"] = [oy, om]
set_classification(r, [("C", "N"), ("C", "P")], ev, "Net new bonds inferred from the 4a-to-5a source drawings; complete mapping is not asserted.")
r["review"]["unresolved"] += ["Purification eluent is not reported for the gram-scale paragraph.", "Candidate structures and partial bond inventory require expert review; SI was not supplied."]
bundle["reactions"].append(r)


# Two separately named alternative phosphorus-source trials, each with only a
# qualitative trace outcome. The paper gives no charge, product structure, or
# measurement method for these trials.
for rid, source, target in [
    ("trace_control_phosphinate", "phosphonate_control", "trace_product_1"),
    ("trace_control_oxaphosphinine", "oxaphosphinine_control", "trace_product_2"),
]:
    ev = [E[3]["text"]]
    r = base_reaction(rid, "specific_experiment", ev, "PDF page 3, paragraph below Scheme 5",
                      "alternative phosphine-source trial", "attempted aminophosphination", ["control"])
    participant(r, "substrate", "4a", "reactant", ev)
    participant(r, "phosphorus_source", source, "reactant", ev)
    participant(r, "trace_product", target, "product", ev)
    s = st.new_step("reported_trial_conditions_unreported", ev)
    s["participant_ids"] = ["substrate", "phosphorus_source", "trace_product"]
    s["operation"] = f("We have tried other phosphine sources", "alternative phosphorus-source trial", ev)
    r["steps"] = [s]
    o = st.new_outcome("trace_product", "trace", [M[target]], ev)
    o["value"] = f("only trace aminophosphination product was detected", "trace", ev)
    o["product_role"] = f("trace aminophosphination product", "trace_product", ev)
    o["scope"] = f("qualitative trace; exact product not identified", {"target": M[target]}, ev)
    r["outcomes"] = [o]
    r["procedure"]["original_text"] = f("We have tried other phosphine sources ... However, only trace aminophosphination product was detected.",
                                         "qualitative control statement", ev)
    set_classification(r, [("C", "N"), ("C", "P")], ev, "Tentative class of the author's stated trace aminophosphination product; its exact structure and mapping are unresolved.")
    r["review"]["unresolved"] += ["Exact product identity, reagent charges, conditions, and trace measurement basis are unreported in pages 1-6.", "No conditions were inherited because the paragraph does not explicitly state that the general procedure was used."]
    bundle["reactions"].append(r)


# Proposed arrow-pushing mechanism is not an experiment.
ev = [E[3]["image"], E[3]["text"]]
r = base_reaction("scheme6_proposed_mechanism", "proposed_mechanism", ev, "PDF page 3, Scheme 6",
                  "Tentative Reaction Mechanism", "proposed electrophilic phosphination/cyclization")
for pid, key, role in [("allene", "1a", "reactant"), ("phosphine_oxide", "2a", "reactant"),
                       ("tf2o", "tf2o", "reagent"), ("lutidine", "lutidine", "reagent"), ("product", "3a", "product")]:
    participant(r, pid, key, role, ev)
r["notes"] = [f("The authors propose electrophilic phosphorus species A, allene-phosphine complex B, intermediates C/D, then deprotonation to 3a.",
                  "author-proposed mechanism with intermediates A-D", ev)]
set_classification(r, [("C", "N"), ("C", "P")], ev, "Net proposed bond formation shown in Scheme 6; this is not experimental proof or a complete elementary-step mapping.")
r["review"]["unresolved"] += ["Intermediates A-D are proposed, not isolated or experimentally verified in pages 1-6.", "DFT details are said to be in unavailable SI."]
bundle["reactions"].append(r)


bundle["source_context"] = [
    {"context_id": f"{DOC}:title_abstract", "document_id": DOC, "evidence_ids": [E[1]["text"], E[1]["image"]],
     "kind": "bibliographic_context", "title": "Metal-Free Intramolecular Aminophosphination of Allenes",
     "doi": "10.1021/acs.joc.0c02169"},
    {"context_id": f"{DOC}:unspecified_initial_catalyst_trials", "document_id": DOC, "evidence_ids": [E[2]["text"]],
     "kind": "unresolved_experiment_group", "description": "Various Ag/Cu catalysts, solvents, and temperatures reportedly gave no phosphinative cyclization product; row-level identities/conditions are not enumerated."},
    {"context_id": f"{DOC}:dfT_discussion", "document_id": DOC, "evidence_ids": [E[3]["text"], E[3]["image"]],
     "kind": "mechanistic_discussion", "description": "DFT/transition-state discussion and author-proposed C-H···π rationale; not a synthetic experiment."},
    {"context_id": f"{DOC}:general_analytical_methods", "document_id": DOC, "evidence_ids": [E[3]["text"]],
     "kind": "analytical_methods", "description": "TLC, NMR, HRMS, chromatography, and allene-substrate preparation citation."},
    {"context_id": f"{DOC}:characterization_pages_4_6", "document_id": DOC,
     "evidence_ids": [E[4]["text"], E[5]["text"], E[6]["text"]], "kind": "characterization",
     "description": "NMR/HRMS characterization blocks for 3a, 3a′, 3b-e, 5a-l, and 6a-d; retained as characterization rather than separate reactions."},
    {"context_id": f"{DOC}:si_mention_only", "document_id": DOC, "evidence_ids": [E[6]["text"], E[6]["image"]],
     "kind": "supplementary_information_notice", "description": "Paper advertises SI, NMR copies, calculation details, and X-ray data for 5g, but no SI file was supplied or associated."},
]

# Raw visual inventory is intentionally separate from the formal bundle.
raw = {
    "batch_id": "allene-stereo-v31-200-20261017", "document_id": DOC,
    "visual_pages_opened_original_detail": [1, 2, 3, 4, 5, 6],
    "paper_title": "Metal-Free Intramolecular Aminophosphination of Allenes",
    "doi": "10.1021/acs.joc.0c02169",
    "decisions": {
        "scheme1": "five prior-literature arrows plus one generic this-work arrow; none counted as a concrete experiment from this paper",
        "scheme2_3": "3a in Scheme 2, Scheme 3, and its product paragraph is one experiment; 3a′/3b-e are separate scope experiments",
        "scheme4": "twelve separate scope experiments; 5k is a hydroamination side-product experiment; 5f/5g dr and 5i/5j E:Z retained",
        "scheme5": "four scope experiments plus two separately named qualitative trace phosphorus-source controls",
        "scheme6": "proposed mechanism, not experimental proof",
        "procedures": "0.2 mmol scope preparations inherit the explicit General Procedures; gram-scale 5a remains separate",
        "characterization": "NMR/HRMS blocks are source_context, not synthetic records",
        "si": "unestablished/not supplied",
    },
    "formal_counts": {"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"])},
}
st.write(RUN / "raw_model/pages_0001_0006_raw_extraction.json", raw)
st.write(RUN / "submitted.json", bundle)

# Stereo diagnostics are persisted for every concrete machine candidate.
diagnostics = {}
for m in bundle["molecules"]:
    smi = st.value(m["smiles"])
    if smi:
        diagnostics[m["molecule_id"]] = st.stereo_inspect(smi)
st.write(RUN / "raw_model/pages_0001_0006_standard_stereo_checks.json", diagnostics)

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "stereo_checks": len(diagnostics)}, ensure_ascii=False))
