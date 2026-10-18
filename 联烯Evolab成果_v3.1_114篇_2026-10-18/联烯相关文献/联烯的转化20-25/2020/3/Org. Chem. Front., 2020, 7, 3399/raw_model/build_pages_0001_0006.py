#!/usr/bin/env python3
"""Build source-grounded v3.1 machine records for PDF pages 1--6.

The script intentionally keeps structures unresolved where the journal drawings
cannot be transcribed safely.  It does not promote any record beyond
``needs_review`` and does not associate the unavailable ESI.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st  # noqa: E402


DOC = "doc_3f18b691f7dce01d"
E = {n: f"{DOC}:page{n}:image" for n in range(1, 7)}
T = {n: f"{DOC}:page{n}:text" for n in range(1, 7)}


def F(raw=None, norm=None, ev=(), method="original_extraction", status=None,
      derivation=None, alternatives=()):
    return st.fact(raw, norm, ev, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def D(norm, ev, why, raw=None):
    return F(raw, norm, ev, method="structure_derivation", derivation=why)


def I(norm, ev, why, raw=None):
    return F(raw, norm, ev, method="inference", derivation=why)


def q(number, unit, raw_text=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if raw_text is not None:
        value["raw_text"] = raw_text
    if qualifier is not None:
        value["qualifier"] = qualifier
    return value


dataset_path = ROOT / "runs/current/dataset.json"
dataset = st.read(dataset_path)
doc = next(d for d in dataset["documents"] if d["document_id"] == DOC)
doc["title"] = ("Phosphine-catalyzed [3 + 2] cycloadditions of trifluoromethyl "
                "enynes/enediynes with allenoates: access to cyclopentenes "
                "containing a CF3-substituted quaternary carbon center")
doc["supplementary_material"] = {
    "status": "unknown",
    "note": ("Page 1 states that ESI is available, but no SI file was supplied "
             "or explicitly associated with this session."),
}
st.write(dataset_path, dataset)

template_path = ROOT / "runs/current/submitted.json"
bundle = st.read(template_path)
bundle["documents"] = copy.deepcopy(dataset["documents"])
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []

molecules = {}


def set_stereo(m, ev, chirality, geometry, representation, *, centers=None,
               double_bonds=None, axial=None, mixture=None, assignment=None,
               isomeric_smiles=None, identity_key=None, why="Source drawing review."):
    s = m["stereochemistry"]
    s["atom_centers"] = D(centers, ev, why, raw="source drawing") if centers is not None else F()
    s["double_bonds"] = D(double_bonds, ev, why, raw="source drawing") if double_bonds is not None else F()
    s["axial"] = D(axial, ev, why, raw="source drawing") if axial is not None else F()
    s["mixture"] = F(mixture, mixture, ev) if mixture is not None else F()
    s["assignment"] = F(assignment, assignment, ev) if assignment is not None else F()
    s["chirality_status"] = I(chirality, ev, why)
    s["geometry_status"] = I(geometry, ev, why)
    s["representation_status"] = I(representation, ev, why)
    s["isomeric_smiles"] = D(isomeric_smiles, ev, why) if isomeric_smiles else F()
    s["identity_key"] = D(identity_key, ev, why) if identity_key else F()


def mol(mid, label, ev, name=None, smiles=None, stereo="unresolved", note=None):
    ev = list(dict.fromkeys(ev))
    if mid in molecules:
        occurrence = molecules[mid]["occurrences"][0]
        occurrence["evidence_ids"] = list(dict.fromkeys(occurrence["evidence_ids"] + ev))
        occurrence["compound_label"]["evidence_ids"] = occurrence["evidence_ids"]
        return mid
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name or label, name or label, ev)
    m["names"]["english"] = F(name, name, ev) if name else F()
    m["names"]["abbreviations"] = F(label, [label], ev) if label != (name or label) else F()
    m["names"]["aliases"] = F()
    m["cas_numbers"] = F()
    m["external_properties"] = []
    m["identity_note"] = F(note, note, ev) if note else F()
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = D(smiles, ev, "Candidate connectivity transcribed from the legible source drawing; RDKit parsing is not source verification.", raw="source structure drawing")
        m["generic_template"] = None
        representation = "partial" if any(p["specified"] == "Unspecified" for p in diag["potential_stereo"]) else "encoded"
        set_stereo(
            m, ev, "achiral", "not_applicable", representation,
            centers=[], double_bonds=[], axial=[],
            isomeric_smiles=diag["canonical_isomeric_smiles"],
            identity_key=diag["stereo_graph_key"],
            why=("The transcribed graph has no potential stereogenic atom, allene axis, "
                 "or configurational E/Z element; terminal allene/alkene termini are CH2 where applicable."),
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F()
        m["generic_template"] = None
        if stereo == "product_unspecified":
            set_stereo(
                m, ev, "unspecified", "not_applicable", "not_encoded",
                centers=[{"feature": "CF3-bearing quaternary cyclopentene carbon", "assignment": "unspecified"}],
                why=("The drawn cyclopentene has a potentially stereogenic CF3-bearing quaternary carbon, "
                     "but no wedge/dash or absolute assignment is printed; no racemic claim is inferred."),
            )
        elif stereo == "diastereomeric_mixture":
            set_stereo(
                m, ev, "mixture", "not_applicable", "not_encoded",
                centers=[{"feature": "new stereocenter(s)", "assignment": "unassigned"}],
                mixture={"type": "diastereomeric mixture", "ratio": "2:1", "configurations": "unassigned"},
                assignment="Relative and absolute configurations are not assigned in the covered pages.",
                why="Scheme 2 reports dr = 2:1 but does not assign the diastereomers.",
            )
        elif stereo == "possible_axial":
            set_stereo(
                m, ev, "unresolved", "unresolved", "unresolved",
                axial=[{"feature": "substituted allene axis", "assignment": "not stated"}],
                why="The substituted allene drawing may support axial stereochemistry, but the covered source gives no assignment.",
            )
        else:
            set_stereo(
                m, ev, "unresolved", "unresolved", "unresolved",
                why="Identity is named or drawn, but a safe complete graph/stereochemical transcription was not established.",
            )
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


def generic_mol(mid, label, ev, fixed, variables, note):
    ev = list(ev)
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(label, label, ev)
    m["names"]["english"] = F()
    m["names"]["abbreviations"] = F()
    m["names"]["aliases"] = F()
    m["cas_numbers"] = F()
    m["representation_type"] = "generic_template"
    m["smiles"] = F()
    m["generic_template"] = {
        "fixed_structure": F(fixed, fixed, ev),
        "variables": [
            {
                "label": v,
                "attachment_points": F(f"source-marked attachment for {v}", f"source-marked attachment for {v}", ev),
                "definition": F("not fully enumerated on the covered page", "not fully enumerated on the covered page", ev),
                "allowed_values": F(),
                "constraints": F(),
            }
            for v in variables
        ],
        "instance_mapping": [],
        "reason": F(note, note, ev),
    }
    set_stereo(m, ev, "unresolved", "unresolved", "not_encoded",
               why="Generic R-group drawing does not identify one concrete stereoisomer.")
    m["identity_note"] = F(note, note, ev)
    m["external_properties"] = []
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


def participant(pid, mid, role, ev, *, amount=None, mass=None, volume=None,
                equiv=None, molpct=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, ev)
    if product_role is not None:
        p["product_role"] = F(product_role, product_role, ev)
    if amount is not None:
        p["charges"]["amount"] = F(amount[0], amount[1], ev)
    if mass is not None:
        p["charges"]["mass"] = F(mass[0], mass[1], ev)
    if volume is not None:
        p["charges"]["volume"] = F(volume[0], volume[1], ev)
    if equiv is not None:
        raw, norm, derived = equiv
        p["charges"]["equivalents"] = I(norm, ev, derived, raw=raw)
    if molpct is not None:
        p["charges"]["mol_percent"] = F(molpct[0], molpct[1], ev)
    if reference is not None:
        p["equivalence_reference"] = F(reference, reference, ev)
    return p


def classification(ev, mode="cycloaddition"):
    ev = list(ev)
    c = st.new_reaction(f"{DOC}:tmp", DOC, ev)["classification"]
    c["mapped_reaction_smiles"] = F()
    c["model_votes"] = []
    if mode == "none":
        c["bond_changes"] = []
        c["labels"] = D([], ev, "No observed transformation; no expected-product bond label is asserted.")
        c["coverage"] = I("not_applicable", ev, "No reaction was observed.")
        c["limitations"] = ["Negative result; intended product connectivity is not classified."]
        return c
    if mode == "unresolved":
        c["bond_changes"] = []
        c["labels"] = D([], ev, "The source reports a complex mixture without an identified product.")
        c["coverage"] = I("unresolved", ev, "Product identities and bond changes are unresolved.")
        c["limitations"] = ["Complex mixture; no product graph or bond mapping is available."]
        return c
    if mode == "c_n":
        changes = [{"change_id": "b1", "atoms": [{"element": "C", "source_role": "alcohol carbon"}, {"element": "N", "source_role": "phthalimide nitrogen"}], "before_order": 0, "after_order": 1, "change_type": "formed"}]
    elif mode == "o_h":
        changes = [
            {"change_id": "b1", "atoms": [{"element": "O", "source_role": "carboxyl oxygen 1"}, {"element": "H", "source_role": "unmapped proton"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
            {"change_id": "b2", "atoms": [{"element": "O", "source_role": "carboxyl oxygen 2"}, {"element": "H", "source_role": "unmapped proton"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        ]
    elif mode == "c_h":
        changes = [{"change_id": "b1", "atoms": [{"element": "C", "source_role": "reduced ester carbon"}, {"element": "H", "source_role": "unmapped hydride/proton"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1}]
    elif mode == "c_c_c_h":
        changes = [
            {"change_id": "b1", "atoms": [{"element": "C", "source_role": "cyclopentene beta carbon"}, {"element": "C", "source_role": "ethyl carbon"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"element": "C", "source_role": "cyclopentene alpha carbon"}, {"element": "H", "source_role": "unmapped proton"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        ]
    elif mode == "c_n_click":
        changes = [{"change_id": "b1", "atoms": [{"element": "C", "source_role": "alkyne carbon"}, {"element": "N", "source_role": "azide nitrogen"}], "before_order": 0, "after_order": 1, "change_type": "formed"}]
    else:
        changes = [
            {"change_id": "b1", "atoms": [{"element": "C", "source_role": "allene-derived terminus"}, {"element": "C", "source_role": "alkene carbon"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"element": "C", "source_role": "allene-derived terminus"}, {"element": "C", "source_role": "alkene carbon"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        ]
    labels = {}
    for ch in changes:
        if ch["change_type"] in ("formed", "hydrogen_gain", "order_increased"):
            lab = st.bond_label(ch["atoms"][0]["element"], ch["atoms"][1]["element"], ch["after_order"])
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(ch["change_id"])
    c["bond_changes"] = [D(ch, ev, "Partial net bond change inferred from the source substrate/product drawings; exact atom mapping is unresolved.") for ch in changes]
    c["labels"] = D(list(labels.values()), ev, "Taxonomy labels derived only from the listed partial net bond changes.")
    c["coverage"] = I("partial", ev, "No complete atom-mapped reaction was asserted.")
    c["limitations"] = [
        "Partial source-drawing classification; exact atom correspondence and all bond-order changes are unresolved.",
        "Unmapped hydrogen/proton sources are not assigned where applicable.",
    ]
    return c


def reaction(rid, ev, record_type, categories, location, *, reaction_name,
             transformation, class_mode="cycloaddition"):
    ev = list(dict.fromkeys(ev))
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories, categories, ev) if categories is not None else F()
    r["location"] = F(location, location, ev)
    r["reaction_name"] = F(reaction_name, reaction_name, ev)
    r["transformation_type"] = I(transformation, ev, "Normalized transformation description derived from the source caption/drawing.")
    r["classification"] = classification(ev, class_mode)
    r["review"] = {
        "status": "needs_review",
        "method": "model",
        "reviewer": None,
        "checked_fields": [],
        "unresolved": [
            "No explicitly associated SI was supplied; procedure, work-up, purification and characterization details may be missing.",
            "Machine extraction and visual transcription require expert review.",
        ],
    }
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


def set_solvent(r, name, ev, volume=None):
    r["solvent_conditions"]["composition"] = F(name, [name], ev)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = F(volume[0], volume[1], ev)


def step(r, sid, ev, pids, *, operation=None, temp=None, time=None, atmosphere=None):
    s = st.new_step(sid, ev)
    s["participant_ids"] = list(pids)
    if operation is not None:
        s["operation"] = F(operation, operation, ev)
    if temp is not None:
        s["temperature"] = F(temp[0], temp[1], ev)
    if time is not None:
        s["time"] = F(time[0], time[1], ev)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, ev)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, ev, raw, norm, *, yield_type=None,
            basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["value"] = F(raw, norm, ev)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, ev)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, ev)
    if scope is not None:
        o["scope"] = F(scope, scope, ev)
    r["outcomes"].append(o)
    return o


# Common named agents/solvents. Their structures are not imported from external
# knowledge; identities and stereochemistry therefore remain unresolved.
agents = {
    "pph3": ("PPh3", "triphenylphosphine"),
    "p4f": ("(4-FC6H4)3P", "tris(4-fluorophenyl)phosphine"),
    "p4ome": ("(4-MeOC6H4)3P", "tris(4-methoxyphenyl)phosphine"),
    "pph2me": ("PPh2Me", "methyldiphenylphosphine"),
    "pbu3": ("PBu3", "tributylphosphine"),
    "dppf": ("DPPF", "1,1'-bis(diphenylphosphino)ferrocene"),
    "binap": ("BINAP", "2,2'-bis(diphenylphosphino)-1,1'-binaphthyl, configuration unspecified"),
    "dppp": ("DPPP", "1,3-bis(diphenylphosphino)propane"),
    "dabco": ("DABCO", "1,4-diazabicyclo[2.2.2]octane"),
    "toluene": ("Toluene", "toluene"), "thf": ("THF", "tetrahydrofuran"),
    "mecn": ("MeCN", "acetonitrile"), "dcm": ("DCM", "dichloromethane"),
    "etoac": ("EtOAc", "ethyl acetate"), "etoh": ("EtOH", "ethanol"),
    "etmgbr": ("EtMgBr", "ethylmagnesium bromide"),
    "cubrsme2": ("CuBr·SMe2", "copper(I) bromide dimethyl sulfide complex"),
    "hmpa": ("HMPA", "hexamethylphosphoramide"), "tmscl": ("TMSCl", "trimethylsilyl chloride"),
    "dibal": ("DIBAL-H", "diisobutylaluminium hydride"),
    "phthalimide": ("Phthalimide", "phthalimide"), "dead": ("DEAD", "diethyl azodicarboxylate"),
    "naoh": ("NaOH", "sodium hydroxide"), "meoh": ("MeOH", "methanol"),
    "phme": ("PhMe", "toluene"), "bnn3": ("BnN3", "benzyl azide"), "cui": ("CuI", "copper(I) iodide"),
}
for key, (label, name) in agents.items():
    mol(f"{DOC}:mol:agent:{key}", label, [E[2], E[3]], name=name,
        note="Named condition component; no external structure or CAS lookup was performed.")

# Concrete allenoate candidates and the safely legible enyne series 2a--2m.
allenoates = {
    "1a": ("ethyl 2,3-butadienoate", "C=C=CC(=O)OCC"),
    "1b": ("benzyl 2,3-butadienoate", "C=C=CC(=O)OCc1ccccc1"),
    "1c": ("cyclohexyl 2,3-butadienoate", "C=C=CC(=O)OC1CCCCC1"),
    "1d": ("tert-butyl 2,3-butadienoate", "C=C=CC(=O)OC(C)(C)C"),
}
for label, (name, smi) in allenoates.items():
    mol(f"{DOC}:mol:{label}", label, [E[2], E[3]], name=name, smiles=smi,
        note="Candidate connectivity transcribed from Table 1/Table 2; terminal allene is not axially chiral.")

enyne_smiles = {
    "2a": "C=C(C#Cc1ccc(C(=O)OC)cc1)C(F)(F)F",
    "2b": "C=C(C#Cc1ccc(C#N)cc1)C(F)(F)F",
    "2c": "C=C(C#Cc1ccc(C(C)=O)cc1)C(F)(F)F",
    "2d": "C=C(C#Cc1ccc(S(C)(=O)=O)cc1)C(F)(F)F",
    "2e": "C=C(C#Cc1ccc(Cl)cc1)C(F)(F)F",
    "2f": "C=C(C#Cc1ccc(Br)cc1)C(F)(F)F",
    "2g": "C=C(C#Cc1ccccc1)C(F)(F)F",
    "2h": "C=C(C#Cc1ccc(C)cc1)C(F)(F)F",
    "2i": "C=C(C#Cc1cccc([N+](=O)[O-])c1)C(F)(F)F",
    "2j": "C=C(C#Cc1cccc(F)c1)C(F)(F)F",
    "2k": "C=C(C#Cc1cccc(C)c1)C(F)(F)F",
    "2l": "C=C(C#Cc1cccc(OC)c1)C(F)(F)F",
    "2m": "C=C(C#Cc1ccccc1C=O)C(F)(F)F",
}
enyne_desc = {
    "2a": "p-CO2Me phenyl", "2b": "p-CN phenyl", "2c": "p-COMe phenyl",
    "2d": "p-SO2Me phenyl", "2e": "p-Cl phenyl", "2f": "p-Br phenyl",
    "2g": "phenyl", "2h": "p-tolyl", "2i": "m-NO2 phenyl", "2j": "m-F phenyl",
    "2k": "m-Me phenyl", "2l": "m-OMe phenyl", "2m": "o-formylphenyl",
    "2n": "Cl/CN-disubstituted phenyl", "2o": "naphthyl", "2p": "indanone-derived aryl",
    "2q": "thienyl", "2r": "benzofuranyl", "2s": "N-tosylcarbazolyl", "2t": "N-tosylindolyl",
}
for label, desc in enyne_desc.items():
    mol(f"{DOC}:mol:{label}", label, [E[3]] + ([E[2]] if label == "2a" else []),
        name=f"2-trifluoromethyl-1,3-enyne ({desc})", smiles=enyne_smiles.get(label),
        note=("Table 2 substituent identity retained. " +
              ("A candidate SMILES was transcribed from the source drawing." if label in enyne_smiles else
               "Exact fused/heteroaryl attachment was not safely converted to SMILES.")))

# Scheme 1 generic and prior-work arrows.
generic_mol(f"{DOC}:mol:g:allenoate", "generic allenoate 1", [E[2]], "C=C=C-C(=O)OR1", ["R1"], "Scheme 1 R-group allenoate.")
generic_mol(f"{DOC}:mol:g:activated_alkene", "generic activated alkene", [E[2]], "R2-C=C-E", ["R2", "E"], "Scheme 1a; E is drawn as carbonyl/phosphonyl/sulfonyl/cyano/nitro activation.")
generic_mol(f"{DOC}:mol:g:lu_products", "generic alpha/gamma cycloadducts", [E[2]], "substituted cyclopentenes", ["R1", "R2", "E"], "Scheme 1a gives two regioisomeric templates.")
generic_mol(f"{DOC}:mol:g:zhang_substrates", "Zhang prior-work substrates", [E[2]], "gamma-substituted allenoate + beta-CF3 enone", ["R1", "R2"], "Prior-work structures in Scheme 1b.")
generic_mol(f"{DOC}:mol:g:zhang_product", "Zhang prior-work cyclopentene", [E[2]], "CF3 cyclopentene", ["R1", "R2"], "Prior-work product in Scheme 1b.")
generic_mol(f"{DOC}:mol:g:loh_substrates", "Loh prior-work substrates", [E[2]], "aryl allenone + beta-CF3 acrylate", ["Ar", "R1"], "Prior-work structures in Scheme 1b.")
generic_mol(f"{DOC}:mol:g:loh_product", "Loh prior-work cyclopentene", [E[2]], "TMS/CF3 cyclopentene", ["Ar", "R1"], "Prior-work product in Scheme 1b.")
generic_mol(f"{DOC}:mol:g:trost_substrates", "Trost prior-work substrates", [E[2]], "TMM precursor + CF3 enyne", ["R"], "Prior-work structures in Scheme 1c.")
generic_mol(f"{DOC}:mol:g:trost_product", "Trost prior-work cyclopentene", [E[2]], "CF3 cyclopentene", ["R"], "Prior-work product in Scheme 1c.")
generic_mol(f"{DOC}:mol:g:current_enyne", "generic CF3-enyne 2", [E[2]], "CH2=C(CF3)-C#C-R2", ["R2"], "This-work generic substrate in Scheme 1d.")
generic_mol(f"{DOC}:mol:g:current_product", "generic alpha-adduct 3", [E[2]], "CF3/alkynyl cyclopentene ester", ["R1", "R2"], "This-work generic product in Scheme 1d.")

scheme_specs = [
    ("scheme1a_lu_literature", f"{DOC}:mol:g:allenoate", f"{DOC}:mol:g:activated_alkene", f"{DOC}:mol:g:lu_products", "Lu's [3 + 2] cycloadditions", "literature example/general prior-work scheme"),
    ("scheme1b_zhang_literature", f"{DOC}:mol:g:zhang_substrates", None, f"{DOC}:mol:g:zhang_product", "Zhang's work", "literature example"),
    ("scheme1b_loh_literature", f"{DOC}:mol:g:loh_substrates", None, f"{DOC}:mol:g:loh_product", "Loh's work", "literature example"),
    ("scheme1c_trost_literature", f"{DOC}:mol:g:trost_substrates", None, f"{DOC}:mol:g:trost_product", "Trost trimethylenemethane cycloaddition", "literature example"),
]
for rid, a, b, p, name, trans in scheme_specs:
    r = reaction(rid, [E[2]], "literature_example", None, {"page": 2, "region": rid}, reaction_name=name, transformation=trans)
    r["participants"] = [participant("input1", a, "reactant", [E[2]])]
    if b:
        r["participants"].append(participant("input2", b, "reactant", [E[2]]))
    r["participants"].append(participant("product", p, "product", [E[2]], product_role="depicted prior-work product/template"))
    r["review"]["unresolved"].append("Scheme is prior work and does not report a current-paper experimental procedure or yield.")

r_generic = reaction("scheme1d_current_generic", [E[2]], "generic_scheme", None,
                     {"page": 2, "region": "Scheme 1d"},
                     reaction_name="Phosphine-catalyzed [3 + 2] cycloaddition of allenoates onto CF3-enynes",
                     transformation="generic [3 + 2] annulation")
r_generic["participants"] = [
    participant("allenoate", f"{DOC}:mol:g:allenoate", "reactant", [E[2]]),
    participant("enyne", f"{DOC}:mol:g:current_enyne", "reactant", [E[2]]),
    participant("pph3", f"{DOC}:mol:agent:pph3", "catalyst", [E[2]]),
    participant("etoac", f"{DOC}:mol:agent:etoac", "solvent", [E[2]]),
    participant("product", f"{DOC}:mol:g:current_product", "product", [E[2]], product_role="alpha-adduct"),
]
set_solvent(r_generic, "EtOAc", [E[2]])
step(r_generic, "reaction", [E[2]], [p["participant_id"] for p in r_generic["participants"]], operation="generic reaction arrow", temp=("rt", "rt"))

# Table 1 product identities.
mol(f"{DOC}:mol:3aa", "3aa", [E[2], E[3]], name="alpha-adduct 3aa", stereo="product_unspecified")
mol(f"{DOC}:mol:4aa", "4aa", [E[2]], name="gamma-adduct 4aa", stereo="product_unspecified")

table1 = [
    (1, "pph3", "toluene", 76, ">20:1"), (2, "pph3", "thf", 87, ">20:1"),
    (3, "pph3", "mecn", 31, ">20:1"), (4, "pph3", "dcm", 66, ">20:1"),
    (5, "pph3", "etoac", 89, ">20:1"), (6, "pph3", "etoh", "Trace", None),
    (7, "p4f", "etoac", 40, ">20:1"), (8, "p4ome", "etoac", 86, ">20:1"),
    (9, "pph2me", "etoac", 84, ">20:1"), (10, "pbu3", "etoac", "Trace", None),
    (11, "dppf", "etoac", 82, ">20:1"), (12, "binap", "etoac", "Trace", None),
    (13, "dppp", "etoac", 57, ">20:1"), (14, "dabco", "etoac", 32, ">20:1"),
]
table1_ids = []
for n, cat, solv, yld, ratio in table1:
    ev = [E[2], E[3]] if n == 5 else [E[2]]
    cats = ["optimization", "substrate_scope"] if n == 5 else ["optimization"]
    r = reaction(f"table1_entry_{n:02d}", ev, "specific_experiment", cats,
                 {"page": 2, "region": f"Table 1 entry {n}"} if n != 5 else
                 {"pages": [2, 3], "regions": ["Table 1 entry 5", "Table 2 product 3aa"]},
                 reaction_name="Phosphine-catalyzed [3 + 2] cycloaddition",
                 transformation="allenoate/enynyl alkene [3 + 2] annulation")
    table1_ids.append(r["reaction_id"])
    r["general_scheme_ids"] = [r_generic["reaction_id"]]
    pids = []
    r["participants"] = [
        participant("allenoate_1a", f"{DOC}:mol:1a", "reactant", ev,
                    amount=("0.4 mmol", q(0.4, "mmol", "0.4 mmol")),
                    equiv=("2.0 equiv (derived from 0.4/0.2 mmol)", q(2.0, "equiv"), "Calculated from explicit 1a and 2a amounts."), reference=f"{DOC}:mol:2a"),
        participant("enyne_2a", f"{DOC}:mol:2a", "reactant", ev,
                    amount=("0.2 mmol", q(0.2, "mmol", "0.2 mmol")),
                    equiv=("1.0 equiv reference", q(1.0, "equiv"), "Chosen denominator is the explicitly smaller 2a charge."), reference=f"{DOC}:mol:2a"),
        participant("catalyst", f"{DOC}:mol:agent:{cat}", "catalyst", ev,
                    molpct=("10 mol%", q(10, "mol%", "10 mol%")), reference=f"{DOC}:mol:2a"),
        participant("solvent", f"{DOC}:mol:agent:{solv}", "solvent", ev,
                    volume=("2 mL", q(2, "mL", "2 mL"))),
        participant("alpha_3aa", f"{DOC}:mol:3aa", "product", ev, product_role="alpha-adduct"),
    ]
    if ratio:
        r["participants"].append(participant("gamma_4aa", f"{DOC}:mol:4aa", "product", ev, product_role="gamma-regioisomer"))
    set_solvent(r, agents[solv][0], ev, ("2 mL", q(2, "mL", "2 mL")))
    r["reaction_scale"] = F("2a (0.2 mmol)", {"reference_molecule_id": f"{DOC}:mol:2a", "amount": q(0.2, "mmol")}, ev)
    step(r, "reaction", ev, [p["participant_id"] for p in r["participants"]],
         operation="Reaction of 1a and 2a with catalyst in solvent",
         temp=("rt", "rt"), time=("24 h", q(24, "h", "24 h")), atmosphere="N2")
    r["procedure"]["original_text"] = F("Conditions: 1a (0.4 mmol), 2a (0.2 mmol), catalyst (10 mol%) in solvent (2 mL) at rt for 24 h under N2.",
                                          "Table 1 general conditions", ev)
    if isinstance(yld, int):
        raw = f"{yld}%" if n != 5 else ["Table 1 entry 5: 89", "Table 2 3aa: 89%"]
        outcome(r, "yield_3aa", "yield", [f"{DOC}:mol:3aa"], ev, raw, q(yld, "%", f"{yld}%"), yield_type="isolated yield", basis="isolated yield", scope="3aa")
    else:
        outcome(r, "trace_3aa", "trace", [f"{DOC}:mol:3aa"], ev, "Trace", "trace", scope="3aa; yield column")
        r["review"]["unresolved"].append("Trace is qualitative; an isolation/measurement method or detection limit is not reported.")
    if ratio:
        outcome(r, "regioselectivity", "selectivity", [f"{DOC}:mol:3aa", f"{DOC}:mol:4aa"], ev,
                ratio, {"metric": "alpha/gamma", "ratio": ratio}, basis="1H NMR analysis of crude products", scope="alpha-adduct 3aa / gamma-adduct 4aa")
    else:
        r["notes"].append(F("Table alpha/gamma cell: —", "table dash (not reported)", ev))
    if n == 10:
        outcome(r, "decomposition", "decomposition", [f"{DOC}:mol:2a"], ev,
                "PBu3 led to rapid decomposition of the starting material", "rapid decomposition of starting material", scope="starting material; extent/method unreported")
    if n == 5:
        r["notes"].append(F("The same 3aa/89% standard-condition experiment appears in Table 1 entry 5 and Table 2.",
                            "duplicate source presentation merged into one experiment", ev, method="inference",
                            derivation="Identical substrates, amounts, catalyst, solvent, temperature, time, atmosphere and yield."))

# Table 2 scope products. Entry 3aa is already Table 1 entry 5.
scope2 = [
    ("3ab", "2b", 70, "16:1", None), ("3ac", "2c", 62, ">20:1", None),
    ("3ad", "2d", 84, ">20:1", None), ("3ae", "2e", 88, ">20:1", None),
    ("3af", "2f", 80, ">20:1", None), ("3ag", "2g", 56, ">20:1", "modified"),
    ("3ah", "2h", 50, ">20:1", "modified"), ("3ai", "2i", 60, ">20:1", None),
    ("3aj", "2j", 65, ">20:1", None), ("3ak", "2k", 73, ">20:1", None),
    ("3al", "2l", 64, ">20:1", None), ("3am", "2m", 75, ">20:1", None),
    ("3an", "2n", 62, ">20:1", None), ("3ao", "2o", 75, ">20:1", None),
    ("3ap", "2p", 51, ">20:1", None), ("3aq", "2q", 70, ">20:1", None),
    ("3ar", "2r", 65, ">20:1", None), ("3as", "2s", 80, ">20:1", None),
    ("3at", "2t", 66, ">20:1", None),
]
table2_ids = [f"{DOC}:table1_entry_05"]
for product, substrate, yld, ratio, modified in scope2:
    pmid = mol(f"{DOC}:mol:{product}", product, [E[3]], name=f"alpha-adduct {product}", stereo="product_unspecified")
    gmid = mol(f"{DOC}:mol:gamma_{product}", f"gamma-{product}", [E[3]], name=f"unlabelled gamma-regioisomer paired with {product}", stereo="product_unspecified")
    r = reaction(f"table2_{product}", [E[3]], "specific_experiment", ["substrate_scope"],
                 {"page": 3, "region": f"Table 2 product {product}"},
                 reaction_name="Phosphine-catalyzed [3 + 2] cycloaddition",
                 transformation="allenoate/enynyl alkene [3 + 2] annulation")
    table2_ids.append(r["reaction_id"])
    r["general_scheme_ids"] = [r_generic["reaction_id"]]
    allene_amt = 0.6 if modified else 0.4
    allene_eq = 3.0 if modified else 2.0
    r["participants"] = [
        participant("allenoate_1a", f"{DOC}:mol:1a", "reactant", [E[3]],
                    amount=(f"{allene_amt} mmol", q(allene_amt, "mmol", f"{allene_amt} mmol")),
                    equiv=(f"{allene_eq} equiv (derived)", q(allene_eq, "equiv"), "Calculated from explicit 1a/2 amounts."), reference=f"{DOC}:mol:{substrate}"),
        participant("enyne", f"{DOC}:mol:{substrate}", "reactant", [E[3]],
                    amount=("0.2 mmol", q(0.2, "mmol", "0.2 mmol")),
                    equiv=("1.0 equiv reference", q(1.0, "equiv"), "Explicit reaction-scale denominator."), reference=f"{DOC}:mol:{substrate}"),
        participant("pph3", f"{DOC}:mol:agent:pph3", "catalyst", [E[3]], molpct=("10 mol%", q(10, "mol%", "10 mol%")), reference=f"{DOC}:mol:{substrate}"),
        participant("etoac", f"{DOC}:mol:agent:etoac", "solvent", [E[3]], volume=("2 mL", q(2, "mL", "2 mL"))),
        participant("alpha_product", pmid, "product", [E[3]], product_role="alpha-adduct"),
        participant("gamma_product", gmid, "product", [E[3]], product_role="unlabelled gamma-regioisomer"),
    ]
    set_solvent(r, "EtOAc", [E[3]], ("2 mL", q(2, "mL", "2 mL")))
    r["reaction_scale"] = F(f"{substrate} (0.2 mmol)", {"reference_molecule_id": f"{DOC}:mol:{substrate}", "amount": q(0.2, "mmol")}, [E[3]])
    time_raw, time_norm = ("48 h", q(48, "h", "48 h")) if modified else ("24 h", q(24, "h", "24 h"))
    step(r, "reaction", [E[3]], [p["participant_id"] for p in r["participants"]], operation="Table 2 scope reaction", temp=("rt", "rt"), time=(time_raw, time_norm), atmosphere="N2")
    r["procedure"]["original_text"] = F("Table 2 general conditions and applicable footnotes", "Table 2 general conditions", [E[3]])
    if modified:
        r["procedure"]["addition_order"] = F("0.6 mmol of 1a was used in two portions", "1a charged in two portions; individual portion timing/size unreported", [E[3]])
        r["procedure"]["modifications"] = F("c 0.6 mmol of 1a was used in two portions; d 48 h", {"allenoate_amount": q(0.6, "mmol"), "time": q(48, "h")}, [E[3]])
    outcome(r, "yield", "yield", [pmid], [E[3]], f"{yld}%", q(yld, "%", f"{yld}%"), yield_type="isolated yield", basis="isolated yield", scope=product)
    outcome(r, "regioselectivity", "selectivity", [pmid, gmid], [E[3]], ratio, {"metric": "alpha/gamma", "ratio": ratio}, scope=f"{product} / unlabelled gamma-regioisomer")
    r["review"]["unresolved"].append("Table 2 does not state the measurement method for the alpha/gamma ratio.")

# Table 2 allenoate variants 3ba--3da.
for product, allene, yld in [("3ba", "1b", 91), ("3ca", "1c", 93), ("3da", "1d", 86)]:
    pmid = mol(f"{DOC}:mol:{product}", product, [E[3]], name=f"alpha-adduct {product}", stereo="product_unspecified")
    gmid = mol(f"{DOC}:mol:gamma_{product}", f"gamma-{product}", [E[3]], name=f"unlabelled gamma-regioisomer paired with {product}", stereo="product_unspecified")
    r = reaction(f"table2_{product}", [E[3]], "specific_experiment", ["substrate_scope"], {"page": 3, "region": f"Table 2 product {product}"}, reaction_name="Phosphine-catalyzed [3 + 2] cycloaddition", transformation="allenoate/enynyl alkene [3 + 2] annulation")
    table2_ids.append(r["reaction_id"])
    r["general_scheme_ids"] = [r_generic["reaction_id"]]
    r["participants"] = [
        participant("allenoate", f"{DOC}:mol:{allene}", "reactant", [E[3]], amount=("0.4 mmol", q(0.4, "mmol")), equiv=("2.0 equiv (derived)", q(2, "equiv"), "Calculated from 0.4/0.2 mmol."), reference=f"{DOC}:mol:2a"),
        participant("enyne_2a", f"{DOC}:mol:2a", "reactant", [E[3]], amount=("0.2 mmol", q(0.2, "mmol")), equiv=("1.0 equiv reference", q(1, "equiv"), "Explicit scale denominator."), reference=f"{DOC}:mol:2a"),
        participant("pph3", f"{DOC}:mol:agent:pph3", "catalyst", [E[3]], molpct=("10 mol%", q(10, "mol%")), reference=f"{DOC}:mol:2a"),
        participant("etoac", f"{DOC}:mol:agent:etoac", "solvent", [E[3]], volume=("2 mL", q(2, "mL"))),
        participant("alpha_product", pmid, "product", [E[3]], product_role="alpha-adduct"),
        participant("gamma_product", gmid, "product", [E[3]], product_role="unlabelled gamma-regioisomer"),
    ]
    set_solvent(r, "EtOAc", [E[3]], ("2 mL", q(2, "mL")))
    r["reaction_scale"] = F("2a (0.2 mmol)", {"reference_molecule_id": f"{DOC}:mol:2a", "amount": q(0.2, "mmol")}, [E[3]])
    step(r, "reaction", [E[3]], [p["participant_id"] for p in r["participants"]], operation="Table 2 allenoate-scope reaction", temp=("rt", "rt"), time=("24 h", q(24, "h")), atmosphere="N2")
    outcome(r, "yield", "yield", [pmid], [E[3]], f"{yld}%", q(yld, "%"), yield_type="isolated yield", basis="isolated yield", scope=product)
    outcome(r, "regioselectivity", "selectivity", [pmid, gmid], [E[3]], ">20:1", {"metric": "alpha/gamma", "ratio": ">20:1"}, scope=f"{product} / gamma regioisomer")
    r["review"]["unresolved"].append("Table 2 does not state the measurement method for the alpha/gamma ratio.")

# Three separately depicted unsuccessful substrates in Table 2.
failed2 = [
    ("alkyl_enyne", "alkyl-substituted CF3-enyne (drawn with -(CH2)2Ph)", "Table 2 unsuccessful substrate, left"),
    ("cf3_styrene", "CF3-styrene", "Table 2 unsuccessful substrate, middle"),
    ("cf3_diene", "CF3-diene", "Table 2 unsuccessful substrate, right"),
]
for key, name, region in failed2:
    mid = mol(f"{DOC}:mol:failed2:{key}", name, [E[3]], name=name,
              note="Unsuccessful substrate drawing retained without inventing a compound number or complete SMILES.")
    r = reaction(f"table2_no_reaction_{key}", [E[3]], "specific_experiment", ["substrate_scope", "control"], {"page": 3, "region": region}, reaction_name="Unsuccessful [3 + 2] cycloaddition substrate", transformation="no observed target cycloaddition", class_mode="none")
    table2_ids.append(r["reaction_id"])
    r["general_scheme_ids"] = [r_generic["reaction_id"]]
    r["participants"] = [
        participant("allenoate_1a", f"{DOC}:mol:1a", "reactant", [E[3]], amount=("0.4 mmol", q(0.4, "mmol")), reference=mid),
        participant("failed_substrate", mid, "reactant", [E[3]], amount=("0.2 mmol", q(0.2, "mmol")), reference=mid),
        participant("pph3", f"{DOC}:mol:agent:pph3", "catalyst", [E[3]], molpct=("10 mol%", q(10, "mol%")), reference=mid),
        participant("etoac", f"{DOC}:mol:agent:etoac", "solvent", [E[3]], volume=("2 mL", q(2, "mL"))),
    ]
    set_solvent(r, "EtOAc", [E[3]], ("2 mL", q(2, "mL")))
    step(r, "reaction", [E[3]], [p["participant_id"] for p in r["participants"]], operation="Table 2 standard-condition unsuccessful-substrate test", temp=("rt", "rt"), time=("24 h", q(24, "h")), atmosphere="N2")
    outcome(r, "no_reaction", "no_reaction", [], [E[3]], "failed to undergo this [3 + 2] cycloaddition with 1a", "no reaction", scope=name)
    r["review"]["unresolved"].append("No detection method, conversion, recovery, or decomposition extent is reported.")

# Table 3 enediyne scope.
mol(f"{DOC}:mol:1e", "1e", [E[3]], name="alpha-trimethylsilyl-substituted aryl allenone 1e", stereo="possible_axial")
for lab, desc in [("1f", "alpha-methyl-substituted allenoate"), ("1g", "alpha-benzyl-substituted allenoate"), ("1h", "gamma-methyl-substituted allenoate"), ("1i", "gamma-phenyl-substituted allenoate")]:
    mol(f"{DOC}:mol:{lab}", lab, [E[3]], name=desc, stereo="possible_axial", note="Substituted allene; axial configuration is not specified.")
for lab, desc in [("5a", "p-bromophenyl CF3-enediyne"), ("5b", "p-tolyl CF3-enediyne"), ("5c", "alkyl/phenethyl CF3-enediyne"), ("5d", "benzyl-ester alkyl CF3-enediyne"), ("5e", "TIPS-ether-containing alkyl CF3-enediyne"), ("5f", "acetone-derived propargylic-alcohol-protected CF3-enediyne"), ("5g", "TIPS-protected CF3-enediyne")]:
    mol(f"{DOC}:mol:{lab}", lab, [E[3]], name=desc, note="Table 3 structure retained; exact complete SMILES and stereochemical graph were not safely transcribed.")

table3_ids = []
for product, sub, allene, yld in [("6a", "5a", "1a", 75), ("6b", "5b", "1a", 92), ("6c", "5c", "1a", 78), ("6d", "5d", "1b", 73), ("6e", "5e", "1a", 77), ("6f", "5f", "1a", 84)]:
    pmid = mol(f"{DOC}:mol:{product}", product, [E[3]], name=f"alpha-adduct {product}", stereo="product_unspecified")
    gmid = mol(f"{DOC}:mol:gamma_{product}", f"gamma-{product}", [E[3]], name=f"unlabelled gamma-regioisomer paired with {product}", stereo="product_unspecified")
    r = reaction(f"table3_{product}", [E[3]], "specific_experiment", ["substrate_scope"], {"page": 3, "region": f"Table 3 product {product}"}, reaction_name="Phosphine-catalyzed [3 + 2] cycloaddition of CF3-enediyne", transformation="allenoate/enediyne [3 + 2] annulation")
    table3_ids.append(r["reaction_id"])
    r["participants"] = [
        participant("allene", f"{DOC}:mol:{allene}", "reactant", [E[3]], amount=("0.2 mmol", q(0.2, "mmol")), equiv=("2.0 equiv (derived)", q(2, "equiv"), "Calculated from 0.2/0.1 mmol."), reference=f"{DOC}:mol:{sub}"),
        participant("enediyne", f"{DOC}:mol:{sub}", "reactant", [E[3]], amount=("0.1 mmol", q(0.1, "mmol")), equiv=("1.0 equiv reference", q(1, "equiv"), "Explicit reaction scale denominator."), reference=f"{DOC}:mol:{sub}"),
        participant("pph3", f"{DOC}:mol:agent:pph3", "catalyst", [E[3]], molpct=("10 mol%", q(10, "mol%")), reference=f"{DOC}:mol:{sub}"),
        participant("etoac", f"{DOC}:mol:agent:etoac", "solvent", [E[3]], volume=("1 mL", q(1, "mL"))),
        participant("alpha_product", pmid, "product", [E[3]], product_role="alpha-adduct"),
        participant("gamma_product", gmid, "product", [E[3]], product_role="unlabelled gamma-regioisomer"),
    ]
    set_solvent(r, "EtOAc", [E[3]], ("1 mL", q(1, "mL")))
    r["reaction_scale"] = F(f"{sub} (0.1 mmol)", {"reference_molecule_id": f"{DOC}:mol:{sub}", "amount": q(0.1, "mmol")}, [E[3]])
    step(r, "reaction", [E[3]], [p["participant_id"] for p in r["participants"]], operation="Table 3 scope reaction", temp=("rt", "rt"), time=("48 h", q(48, "h")), atmosphere="N2")
    outcome(r, "yield", "yield", [pmid], [E[3]], f"{yld}% yield", q(yld, "%"), yield_type="isolated yield", basis="isolated yield", scope=product)
    outcome(r, "regioselectivity", "selectivity", [pmid, gmid], [E[3]], ">20:1", {"metric": "alpha/gamma", "ratio": ">20:1"}, scope=f"{product} / gamma regioisomer")
    r["review"]["unresolved"].append("Table 3 does not state the measurement method for the alpha/gamma ratio.")

r = reaction("table3_6g_no_reaction", [E[3]], "specific_experiment", ["substrate_scope", "control"], {"page": 3, "region": "Table 3 6g, N.R."}, reaction_name="TIPS-protected enediyne no-reaction test", transformation="no observed target cycloaddition", class_mode="none")
table3_ids.append(r["reaction_id"])
r["participants"] = [
    participant("allenoate_1a", f"{DOC}:mol:1a", "reactant", [E[3]], amount=("0.2 mmol", q(0.2, "mmol")), reference=f"{DOC}:mol:5g"),
    participant("enediyne_5g", f"{DOC}:mol:5g", "reactant", [E[3]], amount=("0.1 mmol", q(0.1, "mmol")), reference=f"{DOC}:mol:5g"),
    participant("pph3", f"{DOC}:mol:agent:pph3", "catalyst", [E[3]], molpct=("10 mol%", q(10, "mol%")), reference=f"{DOC}:mol:5g"),
    participant("etoac", f"{DOC}:mol:agent:etoac", "solvent", [E[3]], volume=("1 mL", q(1, "mL"))),
]
set_solvent(r, "EtOAc", [E[3]], ("1 mL", q(1, "mL")))
step(r, "reaction", [E[3]], [p["participant_id"] for p in r["participants"]], operation="Table 3 no-reaction test", temp=("rt", "rt"), time=("48 h", q(48, "h")), atmosphere="N2")
outcome(r, "no_reaction", "no_reaction", [], [E[3]], "6g: N.R.", "no reaction", scope="TIPS-protected enediyne 5g")
r["review"]["unresolved"].append("No product, conversion, recovery, decomposition extent, or detection method is reported.")

# Other allene classes against 5b.
mol(f"{DOC}:mol:6h", "6h", [E[3]], name="alpha-cycloadduct 6h", stereo="product_unspecified")
mol(f"{DOC}:mol:6i", "6i", [E[3]], name="gamma-cycloadduct 6i", stereo="product_unspecified")
r = reaction("table3_1e_to_6h_6i", [E[3]], "specific_experiment", ["substrate_scope"], {"page": 3, "region": "Table 3 other classes of allenes, 1e"}, reaction_name="Cycloaddition of alpha-TMS aryl allenone 1e with enediyne 5b", transformation="allenone/enediyne [3 + 2] annulation")
table3_ids.append(r["reaction_id"])
r["participants"] = [
    participant("allenone_1e", f"{DOC}:mol:1e", "reactant", [E[3]], amount=("0.2 mmol", q(0.2, "mmol")), reference=f"{DOC}:mol:5b"),
    participant("enediyne_5b", f"{DOC}:mol:5b", "reactant", [E[3]], amount=("0.1 mmol", q(0.1, "mmol")), reference=f"{DOC}:mol:5b"),
    participant("pph3", f"{DOC}:mol:agent:pph3", "catalyst", [E[3]], molpct=("10 mol%", q(10, "mol%")), reference=f"{DOC}:mol:5b"),
    participant("etoac", f"{DOC}:mol:agent:etoac", "solvent", [E[3]], volume=("1 mL", q(1, "mL"))),
    participant("alpha_6h", f"{DOC}:mol:6h", "product", [E[3]], product_role="alpha-adduct"),
    participant("gamma_6i", f"{DOC}:mol:6i", "product", [E[3]], product_role="gamma-adduct"),
]
set_solvent(r, "EtOAc", [E[3]], ("1 mL", q(1, "mL")))
step(r, "reaction", [E[3]], [p["participant_id"] for p in r["participants"]], operation="Table 3 other-allene reaction", temp=("rt", "rt"), time=("48 h", q(48, "h")), atmosphere="N2")
outcome(r, "yield_6h", "yield", [f"{DOC}:mol:6h"], [E[3]], "6h: 48% yield", q(48, "%"), yield_type="isolated yield", basis="isolated yield", scope="alpha-adduct 6h")
outcome(r, "yield_6i", "yield", [f"{DOC}:mol:6i"], [E[3]], "6i: 32% yield", q(32, "%"), yield_type="isolated yield", basis="isolated yield", scope="gamma-adduct 6i")
r["notes"].append(F("Table 3 footnote: 'in all cases, the alpha/gamma ratio is >20:1'; depicted yields: 6h 48%, 6i 32%.",
                    "source statements retained without converting isolated yields into a regioisomer ratio", [E[3]]))
r["review"]["unresolved"].append("The blanket Table 3 alpha/gamma >20:1 footnote is inconsistent with the separately depicted 48% alpha and 32% gamma isolated yields; no selected ratio is recorded.")

for lab, result in [("1f", "N.R."), ("1g", "N.R."), ("1h", "complex mixture"), ("1i", "complex mixture")]:
    mode = "none" if result == "N.R." else "unresolved"
    r = reaction(f"table3_{lab}_{'no_reaction' if mode == 'none' else 'complex_mixture'}", [E[3]], "specific_experiment", ["substrate_scope", "control"], {"page": 3, "region": f"Table 3 other classes of allenes, {lab}"}, reaction_name=f"Test of allene {lab} with enediyne 5b", transformation="no target product" if mode == "none" else "unresolved complex mixture", class_mode=mode)
    table3_ids.append(r["reaction_id"])
    r["participants"] = [
        participant("allene", f"{DOC}:mol:{lab}", "reactant", [E[3]], amount=("0.2 mmol", q(0.2, "mmol")), reference=f"{DOC}:mol:5b"),
        participant("enediyne_5b", f"{DOC}:mol:5b", "reactant", [E[3]], amount=("0.1 mmol", q(0.1, "mmol")), reference=f"{DOC}:mol:5b"),
        participant("pph3", f"{DOC}:mol:agent:pph3", "catalyst", [E[3]], molpct=("10 mol%", q(10, "mol%")), reference=f"{DOC}:mol:5b"),
        participant("etoac", f"{DOC}:mol:agent:etoac", "solvent", [E[3]], volume=("1 mL", q(1, "mL"))),
    ]
    set_solvent(r, "EtOAc", [E[3]], ("1 mL", q(1, "mL")))
    step(r, "reaction", [E[3]], [p["participant_id"] for p in r["participants"]], operation="Table 3 unsuccessful/complex-mixture test", temp=("rt", "rt"), time=("48 h", q(48, "h")), atmosphere="N2")
    if mode == "none":
        outcome(r, "no_reaction", "no_reaction", [], [E[3]], "N.R.", "no reaction", scope=lab)
    else:
        outcome(r, "complex_mixture", "other", [], [E[3]], "complex mixture", "complex mixture", scope=lab)
    r["review"]["unresolved"].append("No identified product, measurement method, conversion, or recovery is reported.")

# Scheme 2 scale-up and derivatization arrows.
for lab, name, stereo in [
    ("7", "cyclopentane 7", "diastereomeric_mixture"), ("8", "allylic alcohol 8", "product_unspecified"),
    ("9", "phthalimide-protected allylic amine 9", "product_unspecified"),
    ("10", "dicarboxylic acid 10", "product_unspecified"), ("11", "1,2,3-triazole 11", "product_unspecified"),
]:
    mol(f"{DOC}:mol:{lab}", lab, [E[3], E[4]], name=name, stereo=stereo)

scheme2_ids = []
r = reaction("scheme2_eq1_scale_up_3ae", [E[3], E[4]], "specific_experiment", ["scale_up"], {"pages": [3, 4], "region": "Scheme 2 eqn (1) and continuation prose"}, reaction_name="Gram-scale [3 + 2] cycloaddition", transformation="allenoate/enynyl alkene [3 + 2] annulation")
scheme2_ids.append(r["reaction_id"])
r["general_scheme_ids"] = [r_generic["reaction_id"]]
r["participants"] = [
    participant("allenoate_1a", f"{DOC}:mol:1a", "reactant", [E[3], E[4]], amount=("6 mmol", q(6, "mmol")), mass=("0.67 g", q(0.67, "g")), equiv=("1.5 equivalents", q(1.5, "equiv"), "Source explicitly states 1.5 equivalents."), reference=f"{DOC}:mol:2e"),
    participant("enyne_2e", f"{DOC}:mol:2e", "reactant", [E[3], E[4]], amount=("4 mmol", q(4, "mmol")), mass=("0.92 g", q(0.92, "g")), equiv=("1.0 equiv reference", q(1, "equiv"), "Explicit scale denominator."), reference=f"{DOC}:mol:2e"),
    participant("pph3", f"{DOC}:mol:agent:pph3", "catalyst", [E[3]], molpct=("10 mol%", q(10, "mol%")), reference=f"{DOC}:mol:2e"),
    participant("etoac", f"{DOC}:mol:agent:etoac", "solvent", [E[3]]),
    participant("product_3ae", f"{DOC}:mol:3ae", "product", [E[3], E[4]], product_role="alpha-adduct"),
]
set_solvent(r, "EtOAc", [E[3]])
r["reaction_scale"] = F("2e (0.92 g, 4 mmol)", {"reference_molecule_id": f"{DOC}:mol:2e", "amount": q(4, "mmol"), "mass": q(0.92, "g")}, [E[3], E[4]])
step(r, "reaction", [E[3]], [p["participant_id"] for p in r["participants"]], operation="Gram-scale reaction", temp=("rt", "rt"))
outcome(r, "product_amount", "product_amount", [f"{DOC}:mol:3ae"], [E[3], E[4]], "1.18 g", q(1.18, "g"), scope="3ae")
outcome(r, "yield", "yield", [f"{DOC}:mol:3ae"], [E[3], E[4]], "86%", q(86, "%"), scope="3ae")
r["review"]["unresolved"].append("Scale-up solvent volume, reaction time, atmosphere, work-up, purification, and yield measurement basis are not reported on pages 1–6.")

r = reaction("scheme2_eq2_michael_7", [E[3], E[4]], "specific_experiment", ["derivatization"], {"pages": [3, 4], "region": "Scheme 2 eqn (2)"}, reaction_name="Cu-catalysed Michael addition with ethylmagnesium bromide", transformation="conjugate ethyl addition/protonation", class_mode="c_c_c_h")
scheme2_ids.append(r["reaction_id"])
r["participants"] = [
    participant("substrate_3ae", f"{DOC}:mol:3ae", "reactant", [E[3], E[4]]),
    participant("etmgbr", f"{DOC}:mol:agent:etmgbr", "reagent", [E[3]]), participant("cubrsme2", f"{DOC}:mol:agent:cubrsme2", "catalyst", [E[3]]),
    participant("hmpa", f"{DOC}:mol:agent:hmpa", "additive", [E[3]]), participant("tmscl", f"{DOC}:mol:agent:tmscl", "reagent", [E[3]]),
    participant("thf", f"{DOC}:mol:agent:thf", "solvent", [E[3]]), participant("product_7", f"{DOC}:mol:7", "product", [E[3], E[4]], product_role="cyclopentane product"),
]
set_solvent(r, "THF", [E[3]])
step(r, "reaction", [E[3]], [p["participant_id"] for p in r["participants"]], operation="Michael addition under the displayed reagent set", temp=("-78 °C", q(-78, "°C")))
outcome(r, "yield", "yield", [f"{DOC}:mol:7"], [E[3], E[4]], "62% yield", q(62, "%"), scope="7")
outcome(r, "dr", "selectivity", [f"{DOC}:mol:7"], [E[3], E[4]], "dr = 2:1", {"metric": "dr", "ratio": "2:1", "diastereomer_order": "unassigned"}, scope="diastereomers of 7")
r["review"]["unresolved"].append("Reagent loadings, time, atmosphere, proton source, dr method, and yield measurement basis are unreported.")

r = reaction("scheme2_eq3_reduction_8", [E[3], E[4]], "specific_experiment", ["derivatization"], {"pages": [3, 4], "region": "Scheme 2 eqn (3), first arrow"}, reaction_name="Selective ester reduction with DIBAL-H", transformation="ester-to-alcohol reduction", class_mode="c_h")
scheme2_ids.append(r["reaction_id"])
r["participants"] = [participant("substrate_3ae", f"{DOC}:mol:3ae", "reactant", [E[3], E[4]]), participant("dibal", f"{DOC}:mol:agent:dibal", "reagent", [E[3]]), participant("thf", f"{DOC}:mol:agent:thf", "solvent", [E[3]]), participant("product_8", f"{DOC}:mol:8", "product", [E[3], E[4]], product_role="allylic alcohol")]
set_solvent(r, "THF", [E[3]])
step(r, "reaction", [E[3]], [p["participant_id"] for p in r["participants"]], operation="DIBAL-H reduction")
outcome(r, "yield", "yield", [f"{DOC}:mol:8"], [E[3], E[4]], "89% yield", q(89, "%"), scope="8")
r["review"]["unresolved"].append("Reagent loading, temperature, time, atmosphere, work-up, purification, exact hydrogen mapping, and yield basis are unreported.")

r = reaction("scheme2_eq3_mitsunobu_9", [E[3], E[4]], "specific_experiment", ["derivatization"], {"pages": [3, 4], "region": "Scheme 2 eqn (3), second arrow"}, reaction_name="Mitsunobu reaction with phthalimide", transformation="allylic alcohol-to-protected amine substitution", class_mode="c_n")
scheme2_ids.append(r["reaction_id"])
r["participants"] = [participant("substrate_8", f"{DOC}:mol:8", "reactant", [E[3], E[4]]), participant("phthalimide", f"{DOC}:mol:agent:phthalimide", "reagent", [E[3]]), participant("pph3", f"{DOC}:mol:agent:pph3", "reagent", [E[3]]), participant("dead", f"{DOC}:mol:agent:dead", "reagent", [E[3]]), participant("thf", f"{DOC}:mol:agent:thf", "solvent", [E[3]]), participant("product_9", f"{DOC}:mol:9", "product", [E[3], E[4]], product_role="phthalimide-protected allylic amine")]
set_solvent(r, "THF", [E[3]])
step(r, "reaction", [E[3]], [p["participant_id"] for p in r["participants"]], operation="Mitsunobu substitution", temp=("rt", "rt"))
outcome(r, "yield", "yield", [f"{DOC}:mol:9"], [E[3]], "88% yield", q(88, "%"), scope="9")
r["review"]["unresolved"].append("Reagent loadings, time, atmosphere, work-up, purification, stereochemical course, and yield basis are unreported.")

r = reaction("scheme2_eq4_hydrolysis_10", [E[3], E[4]], "specific_experiment", ["derivatization"], {"pages": [3, 4], "region": "Scheme 2 eqn (4)"}, reaction_name="Hydrolysis to dicarboxylic acid 10", transformation="ester hydrolysis", class_mode="o_h")
scheme2_ids.append(r["reaction_id"])
r["participants"] = [participant("substrate_3ba", f"{DOC}:mol:3ba", "reactant", [E[3], E[4]]), participant("naoh", f"{DOC}:mol:agent:naoh", "reagent", [E[3]]), participant("meoh", f"{DOC}:mol:agent:meoh", "solvent", [E[3]]), participant("product_10", f"{DOC}:mol:10", "product", [E[3], E[4]], product_role="dicarboxylic acid")]
set_solvent(r, "MeOH", [E[3]])
step(r, "reaction", [E[3]], [p["participant_id"] for p in r["participants"]], operation="Hydrolysis with 20% aqueous NaOH", temp=("100 °C", q(100, "°C")))
outcome(r, "yield", "yield", [f"{DOC}:mol:10"], [E[3], E[4]], "95% yield", q(95, "%"), scope="10")
r["notes"].append(F("The structure of compound 10 was unambiguously confirmed by X-ray single-crystal analysis (CCDC 2004999).", {"characterization": "single-crystal X-ray", "CCDC": "2004999"}, [E[4]]))
r["review"]["unresolved"].append("NaOH volume/amount beyond '20% aq.', time, atmosphere, work-up, purification, coproducts, and yield basis are unreported.")

r = reaction("scheme2_eq5_deprotection_click_11", [E[3], E[4]], "specific_experiment", ["derivatization"], {"pages": [3, 4], "region": "Scheme 2 eqn (5), two-step sequence"}, reaction_name="Deprotection followed by azide–alkyne Huisgen cycloaddition", transformation="two-step alkyne deprotection/click cycloaddition", class_mode="c_n_click")
scheme2_ids.append(r["reaction_id"])
r["participants"] = [
    participant("substrate_6f", f"{DOC}:mol:6f", "reactant", [E[3], E[4]]),
    participant("naoh", f"{DOC}:mol:agent:naoh", "reagent", [E[3]], equiv=("2.2 equiv", q(2.2, "equiv"), "Source explicitly reports 2.2 equiv relative to 6f."), reference=f"{DOC}:mol:6f"),
    participant("phme", f"{DOC}:mol:agent:phme", "solvent", [E[3]]),
    participant("bnn3", f"{DOC}:mol:agent:bnn3", "reagent", [E[3]]),
    participant("cui", f"{DOC}:mol:agent:cui", "catalyst", [E[3]], molpct=("50 mol%", q(50, "mol%")), reference=f"{DOC}:mol:6f"),
    participant("meoh", f"{DOC}:mol:agent:meoh", "solvent", [E[3]]),
    participant("product_11", f"{DOC}:mol:11", "product", [E[3], E[4]], product_role="1,2,3-triazole final product"),
]
r["solvent_conditions"]["composition"] = F("i) PhMe; ii) MeOH", [{"stage": "i", "solvent": "PhMe"}, {"stage": "ii", "solvent": "MeOH"}], [E[3]])
step(r, "deprotection", [E[3]], ["substrate_6f", "naoh", "phme"], operation="Treat 6f with NaOH to remove the acetone-derived protecting group", temp=("90 °C", q(90, "°C")), time=("2 h", q(2, "h")))
step(r, "click", [E[3]], ["substrate_6f", "bnn3", "cui", "meoh", "product_11"], operation="Azide–alkyne Huisgen cycloaddition of the resulting terminal alkyne", temp=("80 °C", q(80, "°C")), time=("24 h", q(24, "h")))
outcome(r, "overall_yield", "yield", [f"{DOC}:mol:11"], [E[3], E[4]], "57% yield over two steps", q(57, "%"), yield_type="two-step overall yield", scope={"start": "6f", "end": "11", "steps": 2})
r["procedure"]["addition_order"] = F("i) NaOH (2.2 equiv), PhMe, 90 °C, 2 h; ii) BnN3, CuI (50 mol%), MeOH, 80 °C, 24 h", "stage i then stage ii", [E[3]])
r["review"]["unresolved"].append("Amounts of 6f/BnN3, solvent volumes, atmosphere, work-up, purification, intermediate yield/identity record, and overall-yield measurement basis are unreported.")

# Scheme 3 is an author-proposed mechanism, not a measured sequence.
for lab in ["Int-A", "Int-B", "Int-C", "Int-D", "Int-E"]:
    mol(f"{DOC}:mol:{lab}", lab, [E[4]], name=f"proposed intermediate {lab}", note="Author-proposed mechanistic intermediate; not isolated or experimentally quantified on covered pages.")
r_mech = reaction("scheme3_proposed_mechanism", [E[4]], "proposed_mechanism", ["mechanistic_study"], {"page": 4, "region": "Scheme 3"}, reaction_name="Proposed reaction mechanism", transformation="author-proposed stepwise phosphine catalytic cycle")
r_mech["participants"] = [participant("allenoate_1", f"{DOC}:mol:g:allenoate", "reactant", [E[4]]), participant("enyne_2", f"{DOC}:mol:g:current_enyne", "reactant", [E[4]]), participant("pph3", f"{DOC}:mol:agent:pph3", "catalyst", [E[4]]), participant("product_3", f"{DOC}:mol:g:current_product", "product", [E[4]], product_role="proposed catalytic-cycle product")]
for lab in ["Int-A", "Int-B", "Int-C", "Int-D", "Int-E"]:
    r_mech["participants"].append(participant(lab.lower(), f"{DOC}:mol:{lab}", "unknown", [E[4]]))
r_mech["procedure"]["original_text"] = F("Initially, nucleophilic addition of phosphine onto allenoate 1 ... finally, [1,2] proton transfer followed by beta-elimination ...", "Author-proposed sequence; not an experimental procedure", [E[4]])
r_mech["notes"].append(F("On the basis of the reported literature and our experimental findings, a stepwise mechanism is proposed in Scheme 3.", "author-proposed mechanism", [E[4]]))
r_mech["review"]["unresolved"].append("Mechanistic intermediates and elementary steps are proposed, not directly demonstrated by experiments on pages 1–6.")

# Non-synthetic observations and page inventory context.
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:page1_figure1", "document_id": DOC, "kind": "non_synthetic_figure", "description": "Figure 1 shows selected biologically active cyclopentene/cyclopentane compounds; it is introductory context, not a current-paper reaction experiment.", "evidence_ids": [E[1]]},
    {"context_id": f"{DOC}:context:xray_10", "document_id": DOC, "kind": "characterization", "description": "Compound 10 is reported as confirmed by single-crystal X-ray analysis, CCDC 2004999; no separate synthesis record is created from the characterization statement.", "evidence_ids": [E[4]]},
    {"context_id": f"{DOC}:context:pages5_6_references", "document_id": DOC, "kind": "references", "description": "Pages 5–6 contain bibliography entries only. They were inspected but are not converted into current-paper experiments.", "evidence_ids": [E[5], E[6]]},
]

# Coverage has exactly the requested pages and records every inspected region.
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/abstract and article metadata", "finding": "Current article identity, DOI and summary; page states ESI is available but no SI file was supplied.", "evidence_ids": [E[1], T[1]]},
            {"region": "Figure 1", "finding": "Introductory biological examples; non-synthetic source context, not current experiments.", "evidence_ids": [E[1]]},
        ], "record_ids": [], "unresolved": ["ESI availability is stated, but no SI association/file is available in this session."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Scheme 1a", "finding": "Generic prior-work Lu [3 + 2] cycloaddition; recorded as literature_example.", "evidence_ids": [E[2]]},
            {"region": "Scheme 1b", "finding": "Two distinct prior-work arrows (Zhang and Loh); each recorded separately as literature_example.", "evidence_ids": [E[2]]},
            {"region": "Scheme 1c", "finding": "Trost prior-work TMM cycloaddition; recorded as literature_example.", "evidence_ids": [E[2]]},
            {"region": "Scheme 1d", "finding": "Current-work generic R-group transformation; recorded as generic_scheme.", "evidence_ids": [E[2]]},
            {"region": "Table 1 entries 1–14 and footnotes", "finding": "Fourteen optimization experiments recorded row-by-row; trace/dash outcomes remain qualitative and entry 10 retains the prose decomposition observation.", "evidence_ids": [E[2], T[2]]},
        ], "record_ids": [x["reaction_id"] for x in bundle["reactions"] if x["reaction_id"].endswith(("scheme1a_lu_literature", "scheme1b_zhang_literature", "scheme1b_loh_literature", "scheme1c_trost_literature", "scheme1d_current_generic"))] + table1_ids,
         "unresolved": ["The covered page does not provide work-up/purification details.", "Trace entries lack a detection limit/method; dash selectivity cells remain unreported."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Table 2 products 3aa–3at and 3ba–3da", "finding": "Every depicted scope product mapped; 3aa is the same standard-condition experiment as Table 1 entry 5 and was not duplicated.", "evidence_ids": [E[3], T[3]]},
            {"region": "Table 2 unsuccessful substrates", "finding": "Three separately drawn failed substrate tests recorded without invented products.", "evidence_ids": [E[3], T[3]]},
            {"region": "Table 3 products 6a–6g", "finding": "Six successful scope rows and one no-reaction row recorded.", "evidence_ids": [E[3], T[3]]},
            {"region": "Table 3 other classes of allenes", "finding": "One two-product experiment (6h/6i), two N.R. tests, and two complex-mixture tests recorded separately.", "evidence_ids": [E[3], T[3]]},
            {"region": "Scheme 2 equations (1)–(5)", "finding": "Scale-up plus five derivatization experiment records; eqn (3) split into its two actual arrows, and eqn (5) retained as a two-step overall-yield sequence.", "evidence_ids": [E[3], T[3]]},
        ], "record_ids": list(dict.fromkeys(table2_ids + table3_ids + scheme2_ids)),
         "unresolved": ["Many dense product drawings were not safely converted to complete stereochemical SMILES.", "Table 3 blanket alpha/gamma >20:1 footnote conflicts with depicted 6h 48% and 6i 32% yields; no ratio was selected.", "Scheme 2 omits many reagent doses, times, atmospheres, work-ups, purification details and yield measurement bases."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "continuation of Scheme 2 discussion", "finding": "Text corroborates scale-up masses/yield, derivatization outcomes, dr, two-step overall yield, and X-ray identification of 10.", "evidence_ids": [E[4], T[4]]},
            {"region": "Scheme 3", "finding": "Author-proposed phosphine catalytic cycle; recorded as proposed_mechanism, not as measured experiments.", "evidence_ids": [E[4], T[4]]},
            {"region": "conclusions", "finding": "Summary only; no additional experiment.", "evidence_ids": [E[4], T[4]]},
        ], "record_ids": scheme2_ids + [r_mech["reaction_id"]],
         "unresolved": ["Proposed intermediates are not experimentally quantified on the covered pages."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "references 2–9", "finding": "Bibliography only; cited literature is not imported as current-paper experiments.", "evidence_ids": [E[5], T[5]]},
        ], "record_ids": [], "unresolved": []},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "references 9–16", "finding": "Bibliography only; no current-paper reaction table, arrow or independent characterization dataset.", "evidence_ids": [E[6], T[6]]},
        ], "record_ids": [], "unresolved": []},
    ],
}

# Save raw model extraction separately from the formal submission.  The standard
# merge will quarantine all machine assertions as needs_review.
raw_dir = ROOT / "runs/current/raw_model"
raw_dir.mkdir(parents=True, exist_ok=True)
st.write(raw_dir / "extracted_bundle_pages_0001_0006.json", bundle)
st.write(template_path, bundle)
(ROOT / "runs/current/coverage").mkdir(parents=True, exist_ok=True)
st.write(ROOT / "runs/current/coverage/pages_0001_0006.json", coverage)

# Run and preserve graph-only stereo diagnostics for every concrete candidate.
diagnostics = {}
for m in bundle["molecules"]:
    smi = st.value(m["smiles"])
    if smi:
        diagnostics[m["molecule_id"]] = st.stereo_inspect(smi)
st.write(raw_dir / "stereo_diagnostics.json", diagnostics)

# Literal mention hints are retained only as cross-region review aids; they do
# not merge experiments or elevate verification status.
mentions = st.source_mentions(dataset, ["1a", "2a", "3aa", "3ae", "6f", "N.R.", "complex mixture"], context=220)
st.write(raw_dir / "source_mentions_selected_labels.json", mentions)

print(json.dumps({
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "source_context": len(bundle["source_context"]),
    "concrete_stereo_diagnostics": len(diagnostics),
}, ensure_ascii=False, indent=2))
