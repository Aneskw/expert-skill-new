#!/usr/bin/env python3
"""Build the page-1--9 v3.1 extraction bundle from visually reviewed source pages."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_85bd2ca086151f67"
EV = {p: f"{DOC}:page{p}:image" for p in range(1, 10)}
TXT = {p: f"{DOC}:page{p}:text" for p in range(1, 10)}


def direct(raw, normalized, pages, status=None):
    evidence = [EV[p] for p in pages]
    return st.fact(raw, normalized, evidence, status=status)


def derived(raw, normalized, pages, derivation):
    return st.fact(raw, normalized, [EV[p] for p in pages], method="structure_derivation", derivation=derivation)


def inferred(raw, normalized, pages, derivation):
    return st.fact(raw, normalized, [EV[p] for p in pages], method="inference", derivation=derivation)


def inherited(raw, normalized, pages):
    return st.fact(raw, normalized, [EV[p] for p in pages], method="inherited_general_condition")


def missing(status="unreported"):
    return st.fact(status=status)


def not_applicable():
    return st.fact(method="not_applicable", status="not_applicable")


def set_occurrence(m, label, pages):
    evidence = [EV[p] for p in pages]
    m["occurrences"] = [{
        "document_id": DOC,
        "compound_label": st.fact(label, label, evidence),
        "evidence_ids": evidence,
    }]


def base_molecule(mid, name, label, pages, aliases=()):
    evidence = [EV[p] for p in pages]
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    set_occurrence(m, label, pages)
    m["names"]["original"] = st.fact(name, name, evidence)
    m["names"]["english"] = st.fact(name, name, evidence)
    if aliases:
        m["names"]["aliases"] = st.fact(list(aliases), list(aliases), evidence)
    m["identity_note"] = st.fact(
        "Identity is source-linked but remains a machine extraction pending expert review.",
        "source-linked machine identity; needs review",
        evidence,
    )
    return m


def set_stereo(m, pages, chirality, geometry, representation, *, axial=None, double_bonds=None, mixture=None, assignment=None):
    ev = [EV[p] for p in pages]
    if axial is not None:
        m["stereochemistry"]["axial"] = st.fact(axial, axial, ev)
    if double_bonds is not None:
        m["stereochemistry"]["double_bonds"] = st.fact(double_bonds, double_bonds, ev)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = st.fact(mixture, mixture, ev)
    if assignment is not None:
        m["stereochemistry"]["assignment"] = st.fact(assignment, assignment, ev)
    m["stereochemistry"]["chirality_status"] = st.fact(
        f"Source/structure assessment: {chirality}", chirality, ev,
        method="inference", derivation="Status distinguishes source omission from achirality or racemic composition."
    )
    m["stereochemistry"]["geometry_status"] = st.fact(
        f"Source/structure assessment: {geometry}", geometry, ev,
        method="inference", derivation="Geometry status follows the displayed structure and explicit source wording only."
    )
    m["stereochemistry"]["representation_status"] = st.fact(
        f"Machine representation assessment: {representation}", representation, ev,
        method="inference", derivation="Representation status accounts for unencoded allene/polymer/mixture stereochemistry."
    )


def smiles_molecule(mid, name, label, smiles, pages, *, allene=False, alkene_unspecified=False, aliases=()):
    m = base_molecule(mid, name, label, pages, aliases)
    m["representation_type"] = "smiles"
    m["smiles"] = derived(
        f"Structure/name shown for {label}", smiles, pages,
        "Candidate connectivity was transcribed from the displayed structure/name; RDKit parsing does not verify source identity."
    )
    diag = st.stereo_inspect(smiles)
    if allene:
        geometry = "unspecified" if alkene_unspecified else "not_applicable"
        set_stereo(
            m, pages, "unspecified", geometry, "partial",
            axial={"type": "allene axis", "configuration": None, "source_assignment": "not specified"},
            double_bonds=({"isolated_alkene": "configuration not specified"} if alkene_unspecified else None),
        )
        m["stereochemistry"]["isomeric_smiles"] = derived(
            "Connectivity candidate; allene axial configuration is not encoded", smiles, pages,
            "Toolkit output preserves connectivity but cannot establish unreported allene axial configuration."
        )
    else:
        set_stereo(m, pages, "achiral", "not_applicable", "encoded")
        m["stereochemistry"]["isomeric_smiles"] = derived(
            f"Candidate SMILES for {label}", smiles, pages,
            "No potential tetrahedral or double-bond stereoelement was returned by the required graph diagnostic."
        )
    m["stereochemistry"]["identity_key"] = derived(
        "standard stereo graph diagnostic", diag["stereo_graph_key"], pages,
        "Computed from the candidate isomeric graph; this key is not a sample identity claim."
    )
    return m


def unresolved_molecule(mid, name, label, pages, note=None):
    m = base_molecule(mid, name, label, pages)
    m["representation_type"] = "unresolved"
    set_stereo(m, pages, "unresolved", "unresolved", "unresolved")
    if note:
        m["identity_note"] = direct(note, note, pages)
    return m


def template_fact(raw, normalized, pages, derivation=None):
    if derivation:
        return derived(raw, normalized, pages, derivation)
    return direct(raw, normalized, pages)


def template_molecule(mid, name, label, pages, fixed, *, chirality="unspecified", geometry="not_applicable", representation="not_encoded", mixture=None):
    m = base_molecule(mid, name, label, pages)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": template_fact(
            f"Displayed polymer/template: {fixed}", {"description": fixed}, pages,
            "The source uses a repeat-unit bracket or unresolved macromolecular drawing; no finite molecular SMILES is asserted."
        ),
        "variables": [{
            "label": "n",
            "attachment_points": template_fact("repeat-unit bracket endpoints", ["left", "right"], pages),
            "definition": template_fact("n", "degree of polymerization", pages),
            "allowed_values": missing(),
            "constraints": missing(),
        }],
        "instance_mapping": [],
        "reason": template_fact(
            "Macromolecular repeat unit/end groups are source-defined but chain length and complete molecular identity are not.",
            "polymer represented as a source-linked generic template",
            pages,
        ),
    }
    set_stereo(
        m, pages, chirality, geometry, representation,
        axial={"type": "allene axes in repeat units", "configuration": None, "source_assignment": "not specified"},
        mixture=mixture,
    )
    return m


def set_fact_field(obj, key, raw, normalized, pages, method="direct", derivation=None):
    if method == "inferred":
        obj[key] = inferred(raw, normalized, pages, derivation or "Derived from an explicitly reported ratio.")
    elif method == "inherited":
        obj[key] = inherited(raw, normalized, pages)
    else:
        obj[key] = direct(raw, normalized, pages)


def participant(pid, mid, role, pages, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, [EV[x] for x in pages], stage)
    if product_role is not None:
        p["product_role"] = direct(product_role, product_role, pages)
    return p


def set_charge(p, key, raw, normalized, pages, reference=None, method="direct", derivation=None):
    if method == "inferred":
        p["charges"][key] = inferred(raw, normalized, pages, derivation or "Calculated from the reported feed ratio.")
    elif method == "inherited":
        p["charges"][key] = inherited(raw, normalized, pages)
    else:
        p["charges"][key] = direct(raw, normalized, pages)
    if reference is not None:
        p["equivalence_reference"] = direct(
            f"relative to {reference}", reference, pages
        ) if method != "inferred" else inferred(
            f"ratio denominator is {reference}", reference, pages,
            "The reported feed ratio identifies the loading denominator."
        )


def base_reaction(suffix, pages, record_type="specific_experiment", categories=("other",), reaction_name=None, transformation=None):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, [EV[p] for p in pages])
    r["record_type"] = record_type
    r["experiment_categories"] = direct(list(categories), list(categories), pages)
    r["location"] = direct(
        ", ".join(f"PDF page {p}" for p in pages),
        {"pdf_pages": list(pages)},
        pages,
    )
    if reaction_name:
        r["reaction_name"] = direct(reaction_name, reaction_name, pages)
    if transformation:
        r["transformation_type"] = inferred(
            transformation, transformation, pages,
            "Normalized transformation label inferred from the source reaction name/drawing."
        )
    r["classification"]["coverage"] = inferred(
        "Exact atom mapping not established from the polymer/reaction drawing.", "unresolved", pages,
        "No complete atom mapping is asserted; formation labels remain empty unless explicitly added below."
    )
    r["classification"]["limitations"] = [
        "No automatic atom mapping was performed.",
        "Polymer repeat-unit drawings do not establish a unique finite-molecule mapping.",
    ]
    r["review"]["unresolved"] = [
        "Machine extraction remains needs_review; schema validation is not chemical verification."
    ]
    return r


def add_step(r, sid, pages, *, stage="reaction", pids=(), operation=None, temperature=None, time=None, inherited_condition=False):
    s = st.new_step(sid, [EV[p] for p in pages], stage)
    s["participant_ids"] = list(pids)
    if operation is not None:
        s["operation"] = (inherited(operation, operation, pages) if inherited_condition else direct(operation, operation, pages))
    if temperature is not None:
        raw, norm = temperature
        s["temperature"] = (inherited(raw, norm, pages) if inherited_condition else direct(raw, norm, pages))
    if time is not None:
        raw, norm = time
        s["time"] = (inherited(raw, norm, pages) if inherited_condition else direct(raw, norm, pages))
    r["steps"].append(s)
    return s


def set_solvent(r, name, concentration, pages, *, conflict=None, inherited_condition=False):
    if conflict:
        alternatives = [direct(a[0], a[1], pages) for a in conflict]
        r["solvent_conditions"]["composition"] = st.fact(
            [a[0] for a in conflict], None, [EV[p] for p in pages], status="conflict", alternatives=alternatives
        )
    else:
        r["solvent_conditions"]["composition"] = (
            inherited(name, name, pages) if inherited_condition else direct(name, name, pages)
        )
    if concentration is not None:
        raw, norm = concentration
        r["solvent_conditions"]["concentration"] = (
            inherited(raw, norm, pages) if inherited_condition else direct(raw, norm, pages)
        )
        r["solvent_conditions"]["concentration_basis"] = direct(
            "initial monomer concentration", "initial concentration of cyclic allene monomer", pages
        )


def outcome(r, oid, kind, targets, pages, raw, normalized, *, product_role=None, yield_type=None, basis=None, scope=None, value_fact=None):
    o = st.new_outcome(oid, kind, targets, [EV[p] for p in pages])
    o["value"] = value_fact if value_fact is not None else direct(raw, normalized, pages)
    if product_role is not None:
        o["product_role"] = direct(product_role, product_role, pages)
    if yield_type is not None:
        o["yield_type"] = direct(yield_type, yield_type, pages)
    if basis is not None:
        o["measurement_basis"] = direct(basis, basis, pages)
    if scope is not None:
        o["scope"] = direct(scope, scope, pages)
    r["outcomes"].append(o)
    return o


bundle = json.loads((RUN / "submitted.json").read_text())
fresh_template = json.loads((RUN / "raw_model/init_template_pages_0007_0009.json").read_text())
bundle["documents"] = fresh_template["documents"]
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []

# Molecule registry. CAS values are intentionally unreported because the supplied pages do not give them.
m_cod = smiles_molecule(
    "m_COD", "cis,cis-1,5-cyclooctadiene", "COD", "C1=CCCC=CCC1", [6, 9]
)
set_stereo(
    m_cod, [6, 9], "unspecified", "specified", "partial",
    double_bonds={
        "source_assignment": "cis,cis-1,5-cyclooctadiene",
        "encoding_note": "source geometry is not encoded in the candidate ring SMILES",
    },
)
m_cod["stereochemistry"]["isomeric_smiles"] = derived(
    "Connectivity candidate for COD; the source cis,cis assignment is not encoded",
    "C1=CCCC=CCC1",
    [6, 9],
    "The source names cis,cis-1,5-cyclooctadiene, but the toolkit candidate preserves connectivity only; ring alkene geometry requires expert review.",
)

molecules = [
    smiles_molecule("m_1_cyclononadiene", "1,2-cyclononadiene", "1", "C1=C=CCCCCCC1", [2, 3, 4, 5], allene=True),
    template_molecule("m_poly1", "poly(1)", "poly(1)", [2, 3, 4, 5, 6], "poly(allenamer) repeat unit shown as an opened C9 ring retaining a backbone allene"),
    unresolved_molecule("m_G3", "Grubbs third-generation catalyst", "G3", [2, 3, 4, 5, 6], "Full organoruthenium structure is shown on page 5; no complete coordination-complex SMILES is asserted."),
    unresolved_molecule("m_G2", "Grubbs second-generation catalyst", "G2", [3], "Catalyst label is source-identified; structure is not redrawn on the reviewed page."),
    unresolved_molecule("m_HG2", "Hoveyda–Grubbs second-generation catalyst", "HG2", [3], "Catalyst label is source-identified; structure is not redrawn on the reviewed page."),
    unresolved_molecule("m_G1", "Grubbs first-generation catalyst", "G1", [3], "Catalyst label is source-identified; structure is not redrawn on the reviewed page."),
    unresolved_molecule("m_Mo", "Schrock's catalyst (2,6-diisopropylphenylimidoneophylidene molybdenum(VI) bis(hexafluoro-tert-butoxide))", "Mo", [3], "Catalyst identity is named in prose, but a complete coordination-complex structure is not encoded."),
    smiles_molecule("m_THF", "tetrahydrofuran", "THF", "C1CCOC1", [2, 3, 4]),
    smiles_molecule("m_toluene", "toluene", "toluene", "Cc1ccccc1", [3]),
    smiles_molecule("m_DCM", "dichloromethane", "CH2Cl2", "ClCCl", [3, 6]),
    smiles_molecule("m_hexane", "n-hexane", "n-hexane", "CCCCCC", [3]),
    smiles_molecule("m_chlorobenzene", "chlorobenzene", "chlorobenzene", "Clc1ccccc1", [3]),
    smiles_molecule("m_dioxane", "1,4-dioxane", "1,4-dioxane", "O1CCOCC1", [3]),
    smiles_molecule("m_methanol", "methanol", "methanol", "CO", [3, 4, 6]),
    smiles_molecule("m_EVE", "ethyl vinyl ether", "ethyl vinyl ether", "C=COCC", [3, 4]),
    smiles_molecule("m_CTA2", "2,8-dimethylnona-4,5-diene", "2", "CC(C)CC=C=CCC(C)C", [4], allene=True),
    unresolved_molecule("m_CTA3", "chain-transfer agent 3 (R = CH2Ph in Scheme 2)", "3", [4, 5], "Scheme 2 identifies R = CH2Ph, but the complete finite structure is retained unresolved rather than guessed."),
    unresolved_molecule("m_CTA4", "chain-transfer agent 4 (R = CH2OH in Scheme 2)", "4", [4, 5], "Scheme 2 identifies R = CH2OH, but the complete finite structure is retained unresolved rather than guessed."),
    template_molecule("m_poly1_CTA2", "CTA-2 end-functionalized poly(1)", "poly(1), CTA 2 end groups", [4], "poly(1) bearing isobutyl-derived end groups from CTA 2"),
    template_molecule("m_poly1_CTA3", "CTA-3 end-functionalized poly(1)", "poly(1), CTA 3 end groups", [4, 5], "poly(1) bearing phenethyl end groups from CTA 3"),
    template_molecule("m_poly1_CTA4", "CTA-4 end-functionalized poly(1)", "poly(1), CTA 4 end groups", [4, 5], "poly(1) bearing hydroxyethyl end groups from CTA 4"),
    smiles_molecule("m_CDCl3", "chloroform-d", "CDCl3", "[2H]C(Cl)(Cl)Cl", [5]),
    smiles_molecule(
        "m_THF_d8", "tetrahydrofuran-d8", "THF-d8",
        "[2H]C1([2H])C([2H])([2H])C([2H])([2H])OC1([2H])[2H]", [9]
    ),
    m_cod,
    unresolved_molecule("m_Ru_vinylidene", "ruthenium vinylidene propagating species", "Ru vinylidene", [5], "The organoruthenium polymeric intermediate is shown/assigned, but finite chain length and complete coordination representation are unresolved."),
    smiles_molecule("m_dimethylphenylsilane", "dimethyl(phenyl)silane", "dimethyl(phenyl)silane", "C[SiH](C)c1ccccc1", [6]),
    unresolved_molecule("m_Karstedt", "Karstedt's catalyst", "Karstedt's catalyst", [6], "Named catalyst mixture/complex is source-identified; a unique finite composition is not supplied."),
    template_molecule(
        "m_poly1_Si_mix", "hydrosilylated poly(1) (poly(1)-Si)", "poly(1)-Si", [6],
        "poly(1) hydrosilylation product containing 1,2-addition, E-1,3-addition, and Z-1,3-addition repeat-unit isomers",
        chirality="unspecified", geometry="mixture", representation="not_encoded",
        mixture={"components": ["A: 1,2-addition", "B: E-1,3-addition", "C: Z-1,3-addition"], "source_percentages": [23, 62, 15]},
    ),
    smiles_molecule("m_5_cyclononatriene", "1,2,6-cyclononatriene", "5", "C1=C=CCCC=CCC1", [6], allene=True, alkene_unspecified=True),
    template_molecule(
        "m_poly5", "poly(5)", "poly(5)", [6],
        "polymer with alternating allenyl and alkenyl units from ROMP of 1,2,6-cyclononatriene",
        chirality="unspecified", geometry="mixture", representation="not_encoded",
        mixture={"alkene_geometry": "predominantly cis; exact ratio not reported"},
    ),
]
bundle["molecules"] = molecules

solvent_mid = {
    "THF": "m_THF", "toluene": "m_toluene", "CH2Cl2": "m_DCM", "n-hexane": "m_hexane",
    "chlorobenzene": "m_chlorobenzene", "1,4-dioxane": "m_dioxane",
}
catalyst_mid = {"G3": "m_G3", "G2": "m_G2", "HG2": "m_HG2", "G1": "m_G1", "Mo": "m_Mo"}

# Page 2 current-work generic reaction and Scheme 1 context.
r = base_reaction("scheme1_general_allene_addition_polymerization", [2], "generic_scheme", reaction_name="General allene addition polymerization", transformation="allene addition polymerization")
r["notes"].append(direct("Scheme 1 top uses R = alkyl, alkoxy, amide, aryl, ester derivatives.", "generic R-group allene addition polymerization", [2]))
bundle["reactions"].append(r)

for suffix, name, transformation, note in [
    (
        "scheme1_ochiai_1999_literature",
        "Ochiai (1999) anionic polymerization",
        "literature allene-containing polymer synthesis",
        "Cited Scheme 1 example: anionic polymerization of 4-phenyl-1-buten-3-yne; no present-paper experimental conditions are assigned.",
    ),
    (
        "scheme1_kijima_2002_literature",
        "Kijima (2002) nickel-catalyzed polymerization",
        "literature dehalogenative polycondensation",
        "Cited Scheme 1 example: Ni-catalyzed formation of an alternating allenyl/biphenylene polymer; no present-paper experiment is asserted.",
    ),
    (
        "scheme1_brantley_2020_romp_literature",
        "Brantley (2020) ROMP step",
        "literature norbornene ring-opening metathesis polymerization",
        "First arrow of the cited Brantley sequence; retained separately from the two subsequent transformations.",
    ),
    (
        "scheme1_brantley_2020_dichlorocyclopropanation_literature",
        "Brantley (2020) CHCl3/OH− modification step",
        "literature dichlorocyclopropanation of a polymer alkene",
        "Second arrow of the cited Brantley sequence, labeled CHCl3/OH−; no current-paper experimental conditions are assigned.",
    ),
    (
        "scheme1_brantley_2020_skattebol_literature",
        "Brantley (2020) n-BuLi rearrangement step",
        "literature Doering–Moore–Skattebøl conversion to an allene-containing polymer",
        "Third arrow of the cited Brantley sequence, labeled n-BuLi; not a present-paper experiment.",
    ),
]:
    lr = base_reaction(suffix, [2], "literature_example", reaction_name=name, transformation=transformation)
    lr["notes"].append(direct(note, note, [2]))
    lr["review"]["unresolved"].append("Scheme 1 supplies a literature illustration, not a complete experimental procedure or source-linked molecular registry for this paper.")
    bundle["reactions"].append(lr)

r = base_reaction("scheme1_current_work_romp_1", [2], "generic_scheme", reaction_name="ROMP of 1", transformation="ring-opening allene metathesis polymerization")
r["participants"] = [participant("scheme1_monomer1", "m_1_cyclononadiene", "reactant", [2]), participant("scheme1_poly1", "m_poly1", "product", [2], product_role="poly(allenamer) product")]
r["review"]["unresolved"].append("Scheme 1 is a transformation overview, not a separately quantified experiment.")
bundle["reactions"].append(r)

# Table 1: every optimization row is independent. Entry 1 incorporates the prose experiment rather than duplicating it.
table1 = [
    (1,"G3","THF",100,2.0,25,4.3,7.0,1.6,33,35),
    (2,"G3","toluene",100,2.0,25,None,None,None,0,"trace"),
    (3,"G3","CH2Cl2",100,2.0,25,3.7,3.0,1.8,25,30),
    (4,"G3","n-hexane",100,2.0,25,None,None,None,0,"trace"),
    (5,"G3","chlorobenzene",100,2.0,25,None,None,None,0,"trace"),
    (6,"Mo","THF",100,2.0,25,None,None,None,0,0),
    (7,"G1","THF",100,2.0,25,None,None,None,0,0),
    (8,"G2","THF",100,2.0,25,4.2,8.9,2.0,20,34),
    (9,"HG2","THF",100,2.0,25,4.5,8.6,2.1,25,37),
    (10,"G3","THF",100,2.0,-25,None,None,None,0,0),
    (11,"G3","THF",100,2.0,0,1.3,8.0,2.1,10,11),
    (12,"G3","THF",100,2.0,60,7.3,10.2,2.0,55,60),
    (13,"G2","THF",100,2.0,60,8.6,10.4,2.9,63,70),
    (14,"HG2","THF",100,2.0,60,7.7,11.2,2.6,46,63),
    (15,"G3","1,4-dioxane",100,2.0,80,5.1,10.1,2.5,30,42),
    (16,"G3","1,4-dioxane",100,2.0,100,7.1,12.4,2.5,50,58),
    (17,"G3","THF",50,2.0,60,5.1,7.0,2.6,75,83),
    (18,"G3","THF",25,1.0,60,3.1,6.2,2.0,90,99),
    (19,"G3","THF",25,0.1,60,3.1,1.0,1.6,92,99),
    (20,"G3","THF",25,6.7,60,None,None,None,0,"trace"),
]
for entry, cat, solv, ratio, conc, temp, mn_theo, mn_exp, dispersity, yld, conv in table1:
    pages = [2, 3] if entry == 1 else [3]
    r = base_reaction(f"table1_entry_{entry:02d}", pages, categories=("optimization",), reaction_name="ROMP of 1", transformation="ring-opening allene metathesis polymerization")
    pm = participant(f"t1e{entry}_monomer", "m_1_cyclononadiene", "reactant", pages)
    pc = participant(f"t1e{entry}_catalyst", catalyst_mid[cat], "catalyst", pages)
    set_charge(pc, "mol_percent", f"[M]0/[I]0 = {ratio}", {"number":100/ratio,"unit":"mol%"}, pages, "m_1_cyclononadiene", method="inferred", derivation="Catalyst mol% = 100 divided by the reported monomer-to-initiator feed ratio.")
    ps = participant(f"t1e{entry}_solvent", solvent_mid[solv], "solvent", pages)
    pp = participant(f"t1e{entry}_product", "m_poly1", "product", pages, product_role="poly(1)")
    r["participants"] = [pm, pc, ps, pp]
    if entry == 1:
        pq = participant("t1e1_quench", "m_EVE", "workup_agent", pages, stage="workup")
        set_charge(pq, "equivalents", "excess ethyl vinyl ether", "excess", pages, "m_1_cyclononadiene")
        pmeoh = participant("t1e1_methanol", "m_methanol", "workup_agent", pages, stage="workup")
        set_charge(pmeoh, "volume", "excess methanol", "excess", pages)
        r["participants"] += [pq, pmeoh]
    if entry == 20:
        set_solvent(r, solv, ("[1]0 = 6.7 M", {"number":6.7,"unit":"M"}), pages,
                    conflict=[("Table 1 solvent cell: THF", "THF"), ("footnote g and prose: neat", "neat")])
        r["review"]["unresolved"].append("Table 1 lists THF while footnote g and prose describe neat monomer; neither value is selected.")
    else:
        set_solvent(r, solv, (f"[1]0 = {conc} M", {"number":conc,"unit":"M"}), pages)
    add_step(r, "polymerization", pages, pids=[pm["participant_id"], pc["participant_id"], ps["participant_id"]],
             operation="polymerization mixture held under the Table 1 conditions",
             temperature=(f"{temp} °C", {"number":temp,"unit":"°C"}),
             time=("30 min", {"number":30,"unit":"min"}))
    if entry == 1:
        add_step(r, "quench", [3], stage="workup", pids=[pq["participant_id"]], operation="excess ethyl vinyl ether was added to quench the polymerization")
        add_step(r, "precipitation", [3], stage="workup", pids=[pmeoh["participant_id"], pp["participant_id"]], operation="reaction mixture was poured into excess methanol")
        yfact = st.fact(
            ["Table 1: 33", "prose: yield determined to be 35%"], None, [EV[3]], status="conflict",
            alternatives=[direct("Table 1 entry 1: 33%", {"number":33,"unit":"%"}, [3]), direct("prose: 35%", {"number":35,"unit":"%"}, [3])]
        )
        outcome(r, "yield", "yield", ["m_poly1"], [3], None, None, product_role="poly(1)", value_fact=yfact)
        r["notes"].append(direct("Prose reports 35% yield while Table 1 entry 1 reports 33% yield under the same stated conditions.", "yield conflict retained without choosing a value", [3]))
        r["review"]["unresolved"].append("Yield conflict: prose 35% versus Table 1 entry 1 33%.")
    else:
        outcome(r, "yield", "yield", ["m_poly1"], [3], f"{yld}%", {"number":yld,"unit":"%"}, product_role="poly(1)")
    outcome(r, "conversion", "conversion", ["m_1_cyclononadiene"], [3], str(conv), ({"number":conv,"unit":"%"} if isinstance(conv,(int,float)) else conv), basis="1H NMR analysis of an aliquot removed from the reaction mixture", scope="monomer-to-polymer conversion")
    if mn_theo is not None:
        outcome(r, "mn_theoretical", "other", ["m_poly1"], [3], f"Mn,theo = {mn_theo} kDa", f"{mn_theo} kDa", basis="molecular weight of 1 × [M]0/[I]0 × % conversion", scope="theoretical number-average molecular weight")
    if mn_exp is not None:
        outcome(r, "mn_sec", "other", ["m_poly1"], [3], f"Mn,exp = {mn_exp} kDa", f"{mn_exp} kDa", basis="SEC in THF against polystyrene standards", scope="experimental number-average molecular weight")
    if dispersity is not None:
        outcome(r, "dispersity", "other", ["m_poly1"], [3], f"Đ = {dispersity}", dispersity, basis="SEC", scope="molecular-weight dispersity")
    if entry in (13,14):
        r["notes"].append(direct("Table footnote e: broad, bimodal distribution.", "broad bimodal molecular-weight distribution", [3]))
    if entry in (6,7,10):
        outcome(r, "no_reaction", "no_reaction", [], [3], "0% conversion; no reaction observed", True, scope="polymerization")
    r["review"]["unresolved"] += [
        "Absolute charges, atmosphere, and pressure are not reported in the supplied main-text table.",
        "Yield measurement basis is not explicitly assigned in Table 1; it is not silently labeled isolated or NMR yield.",
        "Exact polymer atom mapping and bond-change classification remain unresolved.",
    ]
    bundle["reactions"].append(r)

# A distinct 1.5 M gelation observation in Table 1 footnote f.
r = base_reaction("table1_footnote_f_1p5M_gelation", [3], categories=("optimization",), reaction_name="ROMP of 1", transformation="ring-opening allene metathesis polymerization")
pm = participant("t1f_monomer", "m_1_cyclononadiene", "reactant", [3])
pc = participant("t1f_catalyst", "m_G3", "catalyst", [3])
set_charge(pc, "mol_percent", "row 18 [M]0/[I]0 = 25", {"number":4,"unit":"mol%"}, [3], "m_1_cyclononadiene", method="inferred")
ps = participant("t1f_solvent", "m_THF", "solvent", [3])
pp = participant("t1f_insoluble", "m_poly1", "product", [3], product_role="insoluble product tentatively similar to poly(1)")
r["participants"] = [pm, pc, ps, pp]
set_solvent(r, "THF", ("[1]0 = 1.5 M", {"number":1.5,"unit":"M"}), [3])
add_step(r, "gelation_test", [3], pids=[x["participant_id"] for x in (pm,pc,ps)], operation="reaction mixture became gelatinous", temperature=("60 °C inherited from row 18", {"number":60,"unit":"°C"}), time=("30 min inherited from Table 1 conditions", {"number":30,"unit":"min"}), inherited_condition=True)
outcome(r, "gelation", "other", ["m_poly1"], [3], "reaction mixture became gelatinous", "gelation", basis="visual observation", scope="bulk reaction behavior")
outcome(r, "ftir_similarity", "other", ["m_poly1"], [3], "FT-IR spectrum of insoluble product was similar to poly(1)", "similar to poly(1)", basis="FT-IR", scope="tentative product characterization")
r["review"]["unresolved"].append("Footnote does not report a yield, conversion, or complete workup for the 1.5 M gelation experiment.")
bundle["reactions"].append(r)

# Scheme 2 generic transformation.
r = base_reaction("scheme2_generic_cta_romp", [4], "generic_scheme", reaction_name="Synthesis of poly(1) in the presence of chain-transfer agents 2–4", transformation="chain-transfer ring-opening allene metathesis polymerization")
r["participants"] = [participant("s2_1", "m_1_cyclononadiene", "reactant", [4]), participant("s2_2", "m_CTA2", "additive", [4]), participant("s2_3", "m_CTA3", "additive", [4]), participant("s2_4", "m_CTA4", "additive", [4])]
r["notes"].append(direct("Scheme 2 defines R = iPr (2), CH2Ph (3), CH2OH (4).", {"2":"iPr","3":"CH2Ph","4":"CH2OH"}, [4]))
r["review"]["unresolved"].append("Generic Scheme 2 is not assigned a single yield or unique CTA loading.")
bundle["reactions"].append(r)

# Table 2 rows; entry 4 incorporates the prose example and its workup/characterization.
table2 = [
    (1,100,50,2.0,9.9,2.5,71), (2,50,50,2.0,7.3,2.8,92), (3,20,50,2.0,6.9,2.8,92),
    (4,10,50,2.0,4.7,2.2,99), (5,10,100,2.0,5.5,2.0,50), (6,10,20,1.0,3.8,2.6,96), (7,10,10,1.0,2.2,2.1,99),
]
for entry, mcta, mcat, conc, mn, dispersity, yld in table2:
    r = base_reaction(f"table2_entry_{entry:02d}", [4], categories=("optimization",), reaction_name="ROMP of 1 in the presence of CTA 2", transformation="chain-transfer ring-opening allene metathesis polymerization")
    pm = participant(f"t2e{entry}_monomer", "m_1_cyclononadiene", "reactant", [4])
    pct = participant(f"t2e{entry}_cta2", "m_CTA2", "additive", [4])
    set_charge(pct, "equivalents", f"[1]0/[2]0 = {mcta}", {"number":1/mcta,"unit":"equiv"}, [4], "m_1_cyclononadiene", method="inferred", derivation="CTA equivalents are the reciprocal of the reported monomer-to-CTA ratio.")
    pc = participant(f"t2e{entry}_g3", "m_G3", "catalyst", [4])
    set_charge(pc, "mol_percent", f"[1]0/[G3]0 = {mcat}", {"number":100/mcat,"unit":"mol%"}, [4], "m_1_cyclononadiene", method="inferred", derivation="Catalyst mol% = 100 divided by the reported monomer-to-G3 ratio.")
    ps = participant(f"t2e{entry}_thf", "m_THF", "solvent", [4])
    pp = participant(f"t2e{entry}_product", "m_poly1_CTA2", "product", [4], product_role="CTA-2 end-functionalized poly(1)")
    r["participants"] = [pm,pct,pc,ps,pp]
    if entry == 4:
        pq = participant("t2e4_quench", "m_EVE", "workup_agent", [4], stage="workup")
        set_charge(pq, "equivalents", "excess ethyl vinyl ether", "excess", [4], "m_1_cyclononadiene")
        pmeoh = participant("t2e4_meoh", "m_methanol", "workup_agent", [4], stage="workup")
        set_charge(pmeoh, "volume", "excess methanol", "excess", [4])
        r["participants"] += [pq,pmeoh]
    set_solvent(r, "THF", (f"[1]0 = {conc} M", {"number":conc,"unit":"M"}), [4])
    add_step(r, "polymerization", [4], pids=[pm["participant_id"],pct["participant_id"],pc["participant_id"],ps["participant_id"]], operation="polymerization in the presence of CTA 2", temperature=("60 °C", {"number":60,"unit":"°C"}), time=("30 min", {"number":30,"unit":"min"}), inherited_condition=(entry != 4))
    if entry == 4:
        add_step(r, "quench", [4], stage="workup", pids=[pq["participant_id"]], operation="excess ethyl vinyl ether was added")
        add_step(r, "precipitation", [4], stage="workup", pids=[pmeoh["participant_id"],pp["participant_id"]], operation="mixture was poured into excess methanol; precipitate was filtered and washed with methanol")
    outcome(r, "yield", "yield", ["m_poly1_CTA2"], [4], f"{yld}%", {"number":yld,"unit":"%"}, product_role="end-functionalized polymer")
    if entry == 4:
        mnfact = st.fact(
            ["Table 2: 4.7 kDa", "prose: 5.0 kDa by SEC"], None, [EV[4]], status="conflict",
            alternatives=[direct("Table 2 entry 4: Mn = 4.7 kDa", "4.7 kDa", [4]), direct("prose: Mn = 5.0 kDa by SEC", "5.0 kDa", [4])]
        )
        outcome(r, "mn_sec", "other", ["m_poly1_CTA2"], [4], None, None, basis="SEC in THF against polystyrene standards", scope="number-average molecular weight", value_fact=mnfact)
        outcome(r, "mn_nmr", "other", ["m_poly1_CTA2"], [4], "Mn = 6.6 kDa from end-group integration", "6.6 kDa", basis="1H NMR end-group integration", scope="number-average molecular weight")
        r["notes"].append(direct("Table 2 reports 4.7 kDa while the matching prose example reports 5.0 kDa by SEC.", "SEC Mn conflict retained without choosing", [4]))
        r["review"]["unresolved"].append("SEC Mn conflict: Table 2 4.7 kDa versus prose 5.0 kDa.")
    else:
        outcome(r, "mn_sec", "other", ["m_poly1_CTA2"], [4], f"Mn = {mn} kDa", f"{mn} kDa", basis="SEC in THF against polystyrene standards", scope="number-average molecular weight")
    outcome(r, "dispersity", "other", ["m_poly1_CTA2"], [4], f"Đ = {dispersity}", dispersity, basis="SEC", scope="molecular-weight dispersity")
    if entry in (6,7):
        r["notes"].append(direct("Footnote c: concentration was reduced to 1.0 M because gelation was observed at 2.0 M.", "2.0 M gelation noted; row run at 1.0 M", [4]))
        r["review"]["unresolved"].append("Footnote c does not state enough detail to assign the prior 2.0 M gelation observation to a separate quantified experiment.")
    r["review"]["unresolved"] += [
        "Absolute charges, atmosphere, and pressure are unreported in the supplied pages.",
        "Yield measurement basis is not explicitly stated in Table 2.",
        "Exact polymer atom mapping and bond-change classification remain unresolved.",
    ]
    bundle["reactions"].append(r)

# CTA 3 and CTA 4 examples, continued onto page 5.
for label, mid, product_mid, loading, mn_nmr, mn_sec, end_group in [
    (3,"m_CTA3","m_poly1_CTA3",2,4.5,4.2,"phenethyl"),
    (4,"m_CTA4","m_poly1_CTA4",1,8.9,7.2,"hydroxyethyl"),
]:
    r = base_reaction(f"scheme2_cta{label}_example", [4,5], categories=("other",), reaction_name=f"ROMP of 1 in the presence of CTA {label}", transformation="chain-transfer ring-opening allene metathesis polymerization")
    pm = participant(f"cta{label}_monomer", "m_1_cyclononadiene", "reactant", [4,5])
    pct = participant(f"cta{label}_agent", mid, "additive", [4,5])
    set_charge(pct, "equivalents", "10:1 monomer-to-CTA molar ratio", {"number":0.1,"unit":"equiv"}, [4], "m_1_cyclononadiene", method="inferred")
    pc = participant(f"cta{label}_g3", "m_G3", "catalyst", [4,5])
    set_charge(pc, "mol_percent", f"{loading} mol % G3", {"number":loading,"unit":"mol%"}, [4], "m_1_cyclononadiene")
    ps = participant(f"cta{label}_thf", "m_THF", "solvent", [4])
    pp = participant(f"cta{label}_product", product_mid, "product", [4,5], product_role=f"poly(1) with {end_group} end groups")
    r["participants"] = [pm,pct,pc,ps,pp]
    set_solvent(r, "THF", ("[1]0 = 2.0 M", {"number":2.0,"unit":"M"}), [4])
    add_step(r, "polymerization", [4], pids=[x["participant_id"] for x in (pm,pct,pc,ps)], operation="ROMP under the Scheme 2 example conditions", temperature=("60 °C inherited from the immediately preceding CTA-2 procedure", {"number":60,"unit":"°C"}), time=("30 min inherited from the immediately preceding CTA-2 procedure", {"number":30,"unit":"min"}), inherited_condition=True)
    outcome(r, "mn_nmr", "other", [product_mid], [5], f"Mn = {mn_nmr} kDa by 1H NMR", f"{mn_nmr} kDa", basis="1H NMR", scope="number-average molecular weight")
    outcome(r, "mn_sec", "other", [product_mid], [5], f"Mn = {mn_sec} kDa by SEC", f"{mn_sec} kDa", basis="SEC", scope="number-average molecular weight")
    r["review"]["unresolved"] += [
        "Temperature and time are inherited from the immediately preceding example rather than repeated in the CTA-3/4 sentence.",
        "No product yield, absolute charge, atmosphere, pressure, or complete workup is reported in the supplied pages.",
    ]
    bundle["reactions"].append(r)

# Figure 3 proposed pathways are a proposal, not a reported experiment.
r = base_reaction("figure3_proposed_mechanism", [5], "proposed_mechanism", reaction_name="Proposed ROMP pathways through ruthenium vinylidene or alkylidene intermediates", transformation="proposed ring-opening allene metathesis mechanism")
r["participants"] = [participant("fig3_monomer", "m_1_cyclononadiene", "reactant", [5]), participant("fig3_g3", "m_G3", "catalyst", [5]), participant("fig3_intermediate", "m_Ru_vinylidene", "product", [5], product_role="proposed propagating intermediate")]
r["notes"] += [
    direct("Figure 3 labels Pathway 1 and Pathway 2 as proposed mechanistic pathways.", ["Pathway 1: Ru vinylidene", "Pathway 2: Ru alkylidene"], [5]),
    direct("The authors conclude that observations are consistent with predominant Pathway 1.", "author interpretation; not direct proof of every elementary step", [5]),
]
r["review"]["unresolved"].append("Elementary bond changes are not classified as observed product formation because Figure 3 is a proposed mechanism.")
bundle["reactions"].append(r)

# Stoichiometric 1H NMR mechanistic experiment.
r = base_reaction("mechanistic_1H_NMR_1_to_G3", [5], categories=("mechanistic_study",), reaction_name="Stoichiometric reaction of 1 with G3 monitored by 1H NMR", transformation="formation of a ruthenium vinylidene propagating species")
pm = participant("mech1_monomer", "m_1_cyclononadiene", "reactant", [5])
pg = participant("mech1_g3", "m_G3", "reactant", [5])
set_charge(pg, "equivalents", "equimolar quantity of G3", {"number":1,"unit":"equiv"}, [5], "m_1_cyclononadiene")
ps = participant("mech1_cdcl3", "m_CDCl3", "solvent", [5])
pp = participant("mech1_product", "m_Ru_vinylidene", "product", [5], product_role="assigned ruthenium vinylidene species with terminal styrenyl groups")
r["participants"] = [pm,pg,ps,pp]
set_solvent(r, "CDCl3", ("[1]0 = 20 mM", {"number":20,"unit":"mM"}), [5])
add_step(r, "stoichiometric_reaction", [5], pids=[x["participant_id"] for x in (pm,pg,ps)], operation="equimolar mixture heated then analyzed by 1H NMR", temperature=("60 °C", {"number":60,"unit":"°C"}), time=("5 min", {"number":5,"unit":"min"}))
outcome(r, "cis_styrenyl", "selectivity", ["m_Ru_vinylidene"], [5], "terminal cis styrenyl moieties ~33%", {"number":33,"unit":"%","qualifier":"approx","raw_text":"~33%"}, basis="1H NMR", scope="terminal cis styrenyl fraction")
outcome(r, "trans_styrenyl", "selectivity", ["m_Ru_vinylidene"], [5], "terminal trans styrenyl moieties ~66%", {"number":66,"unit":"%","qualifier":"approx","raw_text":"~66%"}, basis="1H NMR", scope="terminal trans styrenyl fraction")
outcome(r, "vinylidene_signal", "other", ["m_Ru_vinylidene"], [5], "new broad signal at 5.35 ppm attributed to Ru vinylidene", "5.35 ppm", basis="1H NMR", scope="assigned Ru-vinylidene resonance")
r["review"]["unresolved"] += [
    "Absolute volumes/amounts, atmosphere, and pressure are unreported.",
    "The Ru-vinylidene structure is an assignment; complete finite structure and atom mapping remain unresolved.",
]
bundle["reactions"].append(r)

# 13C NMR mechanistic experiment.
r = base_reaction("mechanistic_13C_NMR_5eq_1_to_G3", [5], categories=("mechanistic_study",), reaction_name="Reaction of G3 with five equivalents of 1 monitored by 13C NMR", transformation="formation of a ruthenium vinylidene species")
pg = participant("mech2_g3", "m_G3", "reactant", [5])
pm = participant("mech2_monomer", "m_1_cyclononadiene", "reactant", [5])
set_charge(pm, "equivalents", "5 molar equivalents of 1", {"number":5,"unit":"equiv"}, [5], "m_G3")
ps = participant("mech2_cdcl3", "m_CDCl3", "solvent", [5])
pp = participant("mech2_product", "m_Ru_vinylidene", "product", [5], product_role="assigned ruthenium vinylidene species")
r["participants"] = [pg,pm,ps,pp]
set_solvent(r, "CDCl3", ("[G3]0 = 360 mM", {"number":360,"unit":"mM"}), [5])
r["solvent_conditions"]["concentration_basis"] = direct("[G3]0", "initial G3 concentration", [5])
add_step(r, "stoichiometric_reaction", [5], pids=[x["participant_id"] for x in (pg,pm,ps)], operation="mixture heated and analyzed by 13C NMR", temperature=("60 °C", {"number":60,"unit":"°C"}))
outcome(r, "vinylidene_13C", "other", ["m_Ru_vinylidene"], [5], "new diagnostic signal at 350 ppm", "350 ppm", basis="13C NMR", scope="assigned Ru-vinylidene carbon resonance")
outcome(r, "monomer_conversion", "conversion", ["m_1_cyclononadiene"], [5], "monomer conversion was quantitative", {"number":100,"unit":"%"}, scope="conversion of 1")
outcome(r, "initiation", "other", ["m_G3"], [5], "initiation was quantitative", "quantitative", scope="G3 initiation")
r["review"]["unresolved"] += [
    "Reaction time, absolute amounts/volume, atmosphere, and pressure are not reported.",
    "The sentence does not separately assign a measurement basis to quantitative initiation/conversion; no method is invented.",
    "Complete organometallic atom mapping and bond classification remain unresolved.",
]
bundle["reactions"].append(r)

# Page 6 hydrosilylation.
r = base_reaction("figure5_poly1_hydrosilylation", [6], categories=("derivatization",), reaction_name="Hydrosilylation of poly(1)", transformation="post-polymerization allene hydrosilylation")
ppoly = participant("hydro_poly1", "m_poly1", "reactant", [6])
psil = participant("hydro_silane", "m_dimethylphenylsilane", "reagent", [6])
set_charge(psil, "equivalents", "1.9 equiv", {"number":1.9,"unit":"equiv"}, [6], "m_poly1")
pcat = participant("hydro_karstedt", "m_Karstedt", "catalyst", [6])
set_charge(pcat, "mol_percent", "5 mol %", {"number":5,"unit":"mol%"}, [6], "m_poly1")
psolv = participant("hydro_dcm", "m_DCM", "solvent", [6])
pprod = participant("hydro_product", "m_poly1_Si_mix", "product", [6], product_role="hydrosilylated polymer isomer mixture")
pmeoh = participant("hydro_meoh", "m_methanol", "workup_agent", [6], stage="workup")
r["participants"] = [ppoly,psil,pcat,psolv,pprod,pmeoh]
set_solvent(r, "CH2Cl2", ("[1]0 = 0.2 M", {"number":0.2,"unit":"M"}), [6])
r["solvent_conditions"]["concentration_basis"] = direct("[1]0 for poly(1)", "initial concentration referenced to poly(1) repeat/allene units", [6])
add_step(r, "hydrosilylation", [6], pids=[x["participant_id"] for x in (ppoly,psil,pcat,psolv)], operation="poly(1) treated with silane and Karstedt's catalyst", temperature=("room temperature", "room temperature"), time=("12 h", {"number":12,"unit":"h"}))
add_step(r, "precipitation", [6], stage="workup", pids=[pmeoh["participant_id"],pprod["participant_id"]], operation="reaction mixture was poured into methanol and the precipitate was collected")
outcome(r, "A_1_2_addition", "selectivity", ["m_poly1_Si_mix"], [6], "1,2-addition product (A), 23%", {"number":23,"unit":"%"}, basis="1H NMR", scope="A = 1,2-addition repeat-unit isomer")
outcome(r, "B_E_1_3_addition", "selectivity", ["m_poly1_Si_mix"], [6], "E-1,3-addition product (B), 62%", {"number":62,"unit":"%"}, basis="1H NMR", scope="B = E-1,3-addition repeat-unit isomer")
outcome(r, "C_Z_1_3_addition", "selectivity", ["m_poly1_Si_mix"], [6], "Z-1,3-addition product (C), 15%", {"number":15,"unit":"%"}, basis="1H NMR", scope="C = Z-1,3-addition repeat-unit isomer")
outcome(r, "product_mn", "other", ["m_poly1_Si_mix"], [6], "poly(1)-Si Mn = 6.0 kDa", "6.0 kDa", basis="SEC", scope="number-average molecular weight")
outcome(r, "starting_mn", "other", ["m_poly1"], [6], "starting poly(1) Mn = 5.3 kDa", "5.3 kDa", basis="SEC", scope="starting-material comparison")
changes = [
    {"change_id":"b1","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"Si"}],"before_order":0,"after_order":1,"change_type":"formed"},
    {"change_id":"b2","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"H"}],"before_order":0,"after_order":1,"change_type":"hydrogen_gain","hydrogen_count":1},
]
labels = [
    {"category_number":12,"code":"C-Si","supporting_change_ids":["b1"]},
    {"category_number":1,"code":"C-H","supporting_change_ids":["b2"]},
]
r["classification"]["bond_changes"] = [inferred("net hydrosilylation bond change", c, [6], "Derived at the repeat-unit level from Si-H addition across the allene; atom identities and all order changes are not fully mapped.") for c in changes]
r["classification"]["labels"] = inferred("C-Si and C-H formation", labels, [6], "Taxonomy labels follow the two asserted net formed bonds only.")
r["classification"]["coverage"] = inferred("repeat-unit net bond formation only", "partial", [6], "Product is a mixture of three repeat-unit isomers; exact atom mapping and accompanying bond-order changes are unresolved.")
r["classification"]["limitations"] = ["C-Si and C-H formation are classified at the net repeat-unit level.", "No complete atom mapping, proton-source mapping, or unique product graph is asserted."]
r["review"]["unresolved"] += [
    "Absolute charges/volumes, atmosphere, pressure, isolated product yield, and purification details are unreported.",
    "The polymer contains a reported A/B/C repeat-unit mixture; no finite isomeric polymer SMILES is asserted.",
]
bundle["reactions"].append(r)

# Page 6 ROMP of compound 5.
r = base_reaction("scheme3_romp_5", [6], categories=("substrate_scope",), reaction_name="ROMP of 5 using G3", transformation="ring-opening allene metathesis polymerization")
pm = participant("s3_monomer5", "m_5_cyclononatriene", "reactant", [6])
pc = participant("s3_g3", "m_G3", "catalyst", [6])
set_charge(pc, "mol_percent", "4 mol % G3", {"number":4,"unit":"mol%"}, [6], "m_5_cyclononatriene")
ps = participant("s3_thf", "m_THF", "solvent", [6])
pp = participant("s3_poly5", "m_poly5", "product", [6], product_role="polymer with alternating allenyl/alkenyl units")
r["participants"] = [pm,pc,ps,pp]
set_solvent(r, "THF", ("[5] = 1.0 M", {"number":1.0,"unit":"M"}), [6])
add_step(r, "polymerization", [6], pids=[x["participant_id"] for x in (pm,pc,ps)], operation="compound 5 subjected to G3 under optimized conditions", temperature=("60 °C", {"number":60,"unit":"°C"}), time=("20 h", {"number":20,"unit":"h"}))
outcome(r, "yield", "yield", ["m_poly5"], [6], "70% yield based on the presumed structure", {"number":70,"unit":"%"}, product_role="poly(5)", basis="based on the presumed polymer structure", scope="polymer obtained")
outcome(r, "mn_sec", "other", ["m_poly5"], [6], "Mn = 2.4 kDa", "2.4 kDa", basis="SEC", scope="number-average molecular weight")
outcome(r, "dispersity", "other", ["m_poly5"], [6], "Đ = 2.2", 2.2, basis="SEC", scope="molecular-weight dispersity")
outcome(r, "alkene_geometry", "selectivity", ["m_poly5"], [6], "alkenes in the polymer backbone were predominately cis", "predominantly cis", basis="13C NMR (27.2 ppm) with FT-IR support (727 cm−1)", scope="backbone alkene geometry; exact cis/trans ratio not reported")
r["review"]["unresolved"] += [
    "Absolute charges/volume, atmosphere, pressure, workup, and purification details are unreported in the supplied main text.",
    "Yield is based on a presumed polymer structure; an isolated/NMR yield label is not added.",
    "Exact cis/trans ratio, allene axial configurations, finite polymer structure, and atom mapping remain unresolved.",
]
bundle["reactions"].append(r)

# Page 9, footnote 85: a distinct 24 h neat-monomer experiment.  It is not
# merged with Table 1 entry 20 because the catalyst ratio and reaction time differ.
r = base_reaction(
    "footnote85_neat_1_G3_24h", [9], categories=("optimization",),
    reaction_name="Neat-monomer reaction of 1 with G3",
    transformation="ring-opening allene metathesis polymerization",
)
pm = participant("fn85_monomer1", "m_1_cyclononadiene", "reactant", [9])
set_charge(pm, "amount", "Approximately 12.3 mmol", {"number":12.3,"unit":"mmol","qualifier":"approx","raw_text":"Approximately 12.3 mmol"}, [9])
set_charge(pm, "volume", "1.5 mL", {"number":1.5,"unit":"mL"}, [9])
pc = participant("fn85_g3", "m_G3", "catalyst", [9])
set_charge(pc, "amount", "0.01 mmol", {"number":0.01,"unit":"mmol"}, [9])
set_charge(pc, "mass", "8.8 mg", {"number":8.8,"unit":"mg"}, [9])
pp = participant("fn85_poly1", "m_poly1", "product", [9], product_role="polymer formed from monomer 1")
pchar = participant("fn85_cdcl3", "m_CDCl3", "solvent", [9], stage="characterization")
r["participants"] = [pm, pc, pp, pchar]
r["reaction_scale"] = direct(
    "Approximately 12.3 mmol (1.5 mL) of 1",
    {"reference_molecule_id":"m_1_cyclononadiene","amount":{"number":12.3,"unit":"mmol","qualifier":"approx","raw_text":"Approximately 12.3 mmol"}},
    [9],
)
r["solvent_conditions"]["composition"] = direct("neat monomer", "neat", [9])
add_step(
    r, "neat_reaction", [9], pids=[pm["participant_id"], pc["participant_id"]],
    operation="Approximately 12.3 mmol (1.5 mL) of 1 was treated with 0.01 mmol (8.8 mg) of G3.",
    time=("After 24 h", {"number":24,"unit":"h"}),
)
add_step(
    r, "aliquot_nmr", [9], stage="characterization", pids=[pchar["participant_id"]],
    operation="an aliquot was removed, diluted with CDCl3, and then analyzed by 1H NMR spectroscopy",
)
outcome(
    r, "monomer_conversion", "conversion", ["m_1_cyclononadiene"], [9],
    "approximately 3% of the monomer converted to polymer",
    {"number":3,"unit":"%","qualifier":"approx","raw_text":"approximately 3%"},
    basis="1H NMR spectroscopy", scope="conversion of monomer 1 to polymer",
)
r["notes"].append(direct(
    "The residual reaction mixture was heterogeneous and comprised of a brown liquid along with a gray solid.",
    {"physical_state":"heterogeneous","components_observed":["brown liquid","gray solid"]}, [9]
))
r["review"]["unresolved"] += [
    "Reaction temperature, atmosphere, pressure, yield, workup, and polymer molecular weight are unreported in footnote 85.",
    "The conversion is approximately 3% by 1H NMR; the source does not report an internal standard or detection limit.",
    "Footnote 85 refers to unavailable Figure S5; no SI evidence was imported.",
    "This 24 h experiment is distinct from Table 1 entry 20 (30 min, [1]0/[G3]0 = 25).",
]
bundle["reactions"].append(r)

# Page 9, footnote 98: the two page-5 NMR experiments were repeated in
# THF-d8.  The source reports only that the results were similar, so exact
# spectral values are not copied from the CDCl3 experiments.
r = base_reaction(
    "mechanistic_1H_NMR_1_to_G3_THF_d8", [5, 9], categories=("mechanistic_study",),
    reaction_name="Stoichiometric reaction of 1 with G3 monitored by 1H NMR in THF-d8",
    transformation="formation of a ruthenium vinylidene propagating species",
)
pm = participant("mech1d8_monomer", "m_1_cyclononadiene", "reactant", [5, 9])
pg = participant("mech1d8_g3", "m_G3", "reactant", [5, 9])
set_charge(pg, "equivalents", "equimolar quantity of G3", {"number":1,"unit":"equiv"}, [5, 9], "m_1_cyclononadiene", method="inherited")
ps = participant("mech1d8_thf", "m_THF_d8", "solvent", [9])
pp = participant("mech1d8_product", "m_Ru_vinylidene", "product", [5, 9], product_role="assigned ruthenium vinylidene species")
r["participants"] = [pm, pg, ps, pp]
r["solvent_conditions"]["composition"] = direct("THF-d8", "THF-d8", [9])
r["solvent_conditions"]["concentration"] = inherited("[1]0 = 20 mM for the referenced 1H NMR experiment", {"number":20,"unit":"mM"}, [5, 9])
r["solvent_conditions"]["concentration_basis"] = inherited("[1]0", "initial concentration of 1", [5, 9])
add_step(
    r, "stoichiometric_reaction", [5, 9], pids=[x["participant_id"] for x in (pm,pg,ps)],
    operation="Performing the NMR experiments in THF-d8 afforded similar results.",
    temperature=("60 °C for the referenced 1H NMR experiment", {"number":60,"unit":"°C"}),
    time=("5 min for the referenced 1H NMR experiment", {"number":5,"unit":"min"}),
    inherited_condition=True,
)
o = outcome(
    r, "similar_1H_result", "other", ["m_Ru_vinylidene"], [9],
    "afforded similar results", "similar to the CDCl3 1H NMR experiment",
    scope="qualitative repetition of the page-5 1H NMR mechanistic experiment in THF-d8",
)
o["measurement_basis"] = inherited("the referenced 1H NMR experiment", "1H NMR", [5, 9])
r["review"]["unresolved"] += [
    "Footnote 98 gives only 'similar results'; exact THF-d8 peak positions and cis/trans fractions are not reported on the supplied page and are not copied from CDCl3.",
    "Concentration, stoichiometry, 60 °C, and 5 min are inherited from the explicitly referenced page-5 1H NMR experiment.",
    "Absolute amounts/volume, atmosphere, pressure, and Figure S9 are unavailable.",
]
bundle["reactions"].append(r)

r = base_reaction(
    "mechanistic_13C_NMR_5eq_1_to_G3_THF_d8", [5, 9], categories=("mechanistic_study",),
    reaction_name="Reaction of G3 with five equivalents of 1 monitored by 13C NMR in THF-d8",
    transformation="formation of a ruthenium vinylidene species",
)
pg = participant("mech2d8_g3", "m_G3", "reactant", [5, 9])
pm = participant("mech2d8_monomer", "m_1_cyclononadiene", "reactant", [5, 9])
set_charge(pm, "equivalents", "5 molar equivalents of 1", {"number":5,"unit":"equiv"}, [5, 9], "m_G3", method="inherited")
ps = participant("mech2d8_thf", "m_THF_d8", "solvent", [9])
pp = participant("mech2d8_product", "m_Ru_vinylidene", "product", [5, 9], product_role="assigned ruthenium vinylidene species")
r["participants"] = [pg, pm, ps, pp]
r["solvent_conditions"]["composition"] = direct("THF-d8", "THF-d8", [9])
r["solvent_conditions"]["concentration"] = inherited("[G3]0 = 360 mM for the referenced 13C NMR experiment", {"number":360,"unit":"mM"}, [5, 9])
r["solvent_conditions"]["concentration_basis"] = inherited("[G3]0", "initial G3 concentration", [5, 9])
add_step(
    r, "stoichiometric_reaction", [5, 9], pids=[x["participant_id"] for x in (pg,pm,ps)],
    operation="Performing the NMR experiments in THF-d8 afforded similar results.",
    temperature=("60 °C for the referenced 13C NMR experiment", {"number":60,"unit":"°C"}),
    inherited_condition=True,
)
o = outcome(
    r, "similar_13C_result", "other", ["m_Ru_vinylidene"], [9],
    "afforded similar results", "similar to the CDCl3 13C NMR experiment",
    scope="qualitative repetition of the page-5 13C NMR mechanistic experiment in THF-d8",
)
o["measurement_basis"] = inherited("the referenced 13C NMR experiment", "13C NMR", [5, 9])
r["review"]["unresolved"] += [
    "Footnote 98 gives only 'similar results'; the exact THF-d8 13C resonance and conversion values are not reported on the supplied page and are not copied from CDCl3.",
    "G3 concentration, five-equivalent loading, and 60 °C are inherited from the explicitly referenced page-5 13C NMR experiment; reaction time remains unreported.",
    "Absolute amounts/volume, atmosphere, pressure, and Figure S9 are unavailable.",
]
bundle["reactions"].append(r)

# Page 9, footnote 112: competition/control experiment between allene 1 and COD.
r = base_reaction(
    "footnote112_competition_1_COD", [9], categories=("control",),
    reaction_name="Competition polymerization of 1 and cis,cis-1,5-cyclooctadiene with G3",
    transformation="selective ring-opening allene metathesis polymerization",
)
pm = participant("fn112_monomer1", "m_1_cyclononadiene", "reactant", [9])
pcod = participant("fn112_cod", "m_COD", "reactant", [9])
set_charge(pcod, "equivalents", "equimolar mixture of 1 and COD", {"number":1,"unit":"equiv"}, [9], "m_1_cyclononadiene", method="inferred", derivation="The source explicitly calls the mixture equimolar and reports equal initial concentrations.")
pc = participant("fn112_g3", "m_G3", "catalyst", [9])
pc["charges"]["mol_percent"] = st.fact(
    "1 mol % G3; loading denominator is not stated",
    None,
    [EV[9]],
    status="unreadable",
)
ps = participant("fn112_thf", "m_THF", "solvent", [9])
pp = participant("fn112_poly1", "m_poly1", "product", [9], product_role="exclusive polymer product")
r["participants"] = [pm, pcod, pc, ps, pp]
r["solvent_conditions"]["composition"] = direct("THF", "THF", [9])
r["solvent_conditions"]["concentration"] = direct(
    "[COD]0 = [1]0 = 2.0 M",
    {"m_COD":{"number":2.0,"unit":"M"},"m_1_cyclononadiene":{"number":2.0,"unit":"M"}}, [9]
)
r["solvent_conditions"]["concentration_basis"] = direct(
    "[COD]0 = [1]0", ["initial COD concentration","initial concentration of 1"], [9]
)
add_step(
    r, "competition_polymerization", [9], pids=[x["participant_id"] for x in (pm,pcod,pc,ps)],
    operation="treatment of an equimolar mixture of 1 and COD with G3",
    temperature=("25 °C", {"number":25,"unit":"°C"}),
)
outcome(
    r, "exclusive_poly1", "selectivity", ["m_poly1"], [9],
    "resulted in the formation of poly(1) as the exclusive product", "poly(1) exclusive",
    product_role="exclusive polymer product",
    scope="chemoselectivity between reaction at the allene of 1 and alkene units of COD",
)
r["review"]["unresolved"] += [
    "Reaction time, absolute charges/volume, atmosphere, pressure, conversion, yield, workup, and measurement basis for 'exclusive product' are unreported.",
    "The denominator for the reported 1 mol % G3 is not explicitly assigned, so no equivalence reference is invented.",
    "No COD-derived product is reported, but a quantitative no-reaction outcome for COD is not asserted without a stated method or detection limit.",
    "Exact polymer atom mapping and bond-formation classification remain unresolved.",
]
bundle["reactions"].append(r)

# Source context preserves non-synthetic observations and classifies cited work without promoting it to current experiments.
bundle["source_context"] = [
    {
        "context_id": "ctx_page1_toc_intro",
        "document_id": DOC,
        "evidence_ids": [EV[1]],
        "kind": "article_context",
        "summary": "Title/abstract/graphical abstract and introduction. The graphic summarizes this work; prose cites prior allene polymerizations and transformations, not present-paper experiments on page 1.",
    },
    {
        "context_id": "ctx_page2_scheme1_figure1",
        "document_id": DOC,
        "evidence_ids": [EV[2]],
        "kind": "scheme_and_characterization_context",
        "scheme1_arrow_classification": [
            {"label":"General Allene Addition Polymerization","classification":"generic_scheme"},
            {"label":"Ochiai (1999)","classification":"literature_example"},
            {"label":"Kijima (2002)","classification":"literature_example"},
            {"label":"Brantley (2020)","classification":"literature_example"},
            {"label":"Neary & Moore / This Work (2021)","classification":"generic current-work overview plus cited precedent"},
        ],
        "summary": "Figure 1 is a 1H NMR comparison of monomer 1, poly(1), and the cited addition polymer; it is characterization rather than a separate synthesis record.",
    },
    {
        "context_id": "ctx_page3_poly1_characterization",
        "document_id": DOC,
        "evidence_ids": [EV[3]],
        "kind": "characterization_context",
        "summary": "Poly(1) characterization: 1H NMR 5.04 ppm, 13C NMR 204 ppm, FT-IR 1960 cm−1; monomer/addition-polymer comparison. These observations support Table 1 products and are not separate synthetic reactions.",
    },
    {
        "context_id": "ctx_page4_figure2_table2",
        "document_id": DOC,
        "evidence_ids": [EV[4]],
        "kind": "characterization_and_footnote_context",
        "summary": "Figure 2 is 1H NMR characterization of Table 2 entry 4 material. Footnote c reports gelation at 2.0 M but does not unambiguously define separate quantified runs for entries 6 and 7.",
    },
    {
        "context_id": "ctx_page5_mechanistic_spectra",
        "document_id": DOC,
        "evidence_ids": [EV[5]],
        "kind": "mechanistic_context",
        "summary": "Figure 3 is explicitly proposed mechanism; Figure 4 is the spectrum from the separately recorded stoichiometric 1H NMR mechanistic experiment. Spectral peaks are not standalone syntheses.",
    },
    {
        "context_id": "ctx_page6_characterization_and_si_notice",
        "document_id": DOC,
        "evidence_ids": [EV[6]],
        "kind": "characterization_and_document_context",
        "summary": "Figure 5 NMR/SEC and poly(5) NMR/FT-IR observations are linked to their reaction records. The page says SI exists online, but no SI file was explicitly supplied or associated, so no SI evidence is imported.",
    },
    {
        "context_id": "ctx_page7_author_references",
        "document_id": DOC,
        "evidence_ids": [EV[7]],
        "kind": "document_and_literature_context",
        "summary": "The page finishes the Associated Content description, gives author/acknowledgment information, and lists references 1–31. The bibliography describes cited literature and is not promoted to present-paper experiments; no SI file was supplied.",
    },
    {
        "context_id": "ctx_page8_references",
        "document_id": DOC,
        "evidence_ids": [EV[8]],
        "kind": "literature_context",
        "summary": "The full page contains references 32–75 only. These are citations, with no current-paper reaction arrow, table row, procedure, result, or footnote experiment on this page.",
    },
    {
        "context_id": "ctx_page9_references_and_footnotes",
        "document_id": DOC,
        "evidence_ids": [EV[9]],
        "kind": "literature_characterization_and_footnote_context",
        "summary": "The page contains references 75–112 plus present-paper footnotes. Footnotes 85, 98, and 112 map to specific experiment records. Footnote 93 reports styrenyl end-group signals in isolated poly(1), and footnote 99 points to a full 1H NMR spectrum; these are characterization pointers rather than new syntheses. Unavailable SI figures S5, S7, S9, and S10 were not imported.",
    },
]

# Raw extraction snapshot keeps the visually transcribed rows and page decisions separate from the normalized bundle.
raw = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "reviewed_pdf_pages": [1,2,3,4,5,6],
    "source_title": "Ring Opening Metathesis Polymerization of Cyclic Allenes",
    "doi": "10.1021/acs.macromol.1c00571",
    "table1_rows": table1,
    "table2_rows": table2,
    "visual_conflicts": [
        "Table 1 entry 1 yield 33% versus matching prose yield 35%.",
        "Table 2 entry 4 SEC Mn 4.7 kDa versus matching prose SEC Mn 5.0 kDa.",
        "Table 1 entry 20 solvent cell THF versus footnote/prose neat.",
    ],
    "record_count": 43,
    "molecule_count": 28,
    "semantic_accuracy_verified": False,
}
st.write(RUN / "raw_model/extracted_pages_0001_0006.json", raw)
raw_pages_7_9 = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "reviewed_pdf_pages": [7, 8, 9],
    "source_title": "Ring Opening Metathesis Polymerization of Cyclic Allenes",
    "doi": "10.1021/acs.macromol.1c00571",
    "page_findings": {
        "7": "Associated Content tail, author information, acknowledgments, and references 1-31; no unique present-paper experiment.",
        "8": "References 32-75 only; no present-paper experiment.",
        "9": "References 75-112 plus footnote experiments 85, 98, and 112; footnotes 93 and 99 are characterization/spectrum pointers.",
    },
    "new_reaction_ids": [
        f"{DOC}:footnote85_neat_1_G3_24h",
        f"{DOC}:mechanistic_1H_NMR_1_to_G3_THF_d8",
        f"{DOC}:mechanistic_13C_NMR_5eq_1_to_G3_THF_d8",
        f"{DOC}:footnote112_competition_1_COD",
    ],
    "new_molecule_ids": ["m_THF_d8", "m_COD"],
    "unavailable_references": ["Figure S5", "Figure S7", "Figure S9", "Figure S10"],
    "semantic_accuracy_verified": False,
}
st.write(RUN / "raw_model/extracted_pages_0007_0009.json", raw_pages_7_9)
st.write(RUN / "submitted.json", bundle)
print(json.dumps({"molecules":len(bundle["molecules"]),"reactions":len(bundle["reactions"]),"source_context":len(bundle["source_context"])}, ensure_ascii=False))
