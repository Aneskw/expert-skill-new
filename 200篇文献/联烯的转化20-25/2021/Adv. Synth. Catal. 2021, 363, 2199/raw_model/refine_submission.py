#!/usr/bin/env python3
"""Apply source-backed structure refinements found during advisory audit."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_submission as b


def main():
    bundle = b.st.read(b.SUBMITTED)
    replacements = {
        f"{b.DOC}:mol:1c": b.molecule(
            f"{b.DOC}:mol:1c", "1c", [b.ev(3, "image"), b.ev(3)],
            "yne-allenone with R1 = 4-chlorophenyl",
            "C=C=CC(=O)c1ccccc1C#Cc1ccc(Cl)cc1", allene=True,
            identity_note="Candidate SMILES encoded because both the drawing and prose specify the para-chlorophenyl substituent."),
        f"{b.DOC}:mol:2f": b.molecule(
            f"{b.DOC}:mol:2f", "2f", [b.ev(3, "image"), b.ev(3)],
            "1-alkynylnaphthalen-2-ol with R3 = 4-tert-butylphenyl",
            "Oc1ccc2ccccc2c1C#Cc1ccc(C(C)(C)C)cc1",
            identity_note="Candidate SMILES encoded because both the drawing and prose specify the para-tert-butylphenyl substituent."),
    }
    bundle["molecules"] = [replacements.get(m["molecule_id"], m) for m in bundle["molecules"]]
    scope3a = next(r for r in bundle["reactions"] if r["reaction_id"] == f"{b.DOC}:scheme2_scope_3a")
    direct = {p["participant_id"]: p for p in scope3a["participants"]}
    direct["scope_1"]["charges"]["mass"] = b.fact("48.8 mg", b.q(48.8, "mg"), [b.ev(5)])
    direct["scope_1"]["charges"]["amount"] = b.fact("0.2 mmol", b.q(0.2, "mmol"), [b.ev(5)])
    direct["scope_1"]["charges"]["equivalents"] = b.fact("1 equiv.", b.q(1, "equiv"), [b.ev(5)])
    direct["scope_1"]["equivalence_reference"] = b.fact("1a", f"{b.DOC}:mol:1a", [b.ev(5)])
    direct["scope_2"]["charges"]["mass"] = b.fact("107.36 mg", b.q(107.36, "mg"), [b.ev(5)])
    direct["scope_2"]["charges"]["amount"] = b.fact("0.44 mmol", b.q(0.44, "mmol"), [b.ev(5)])
    direct["scope_2"]["charges"]["equivalents"] = b.fact("2.2 equiv.", b.q(2.2, "equiv"), [b.ev(5)])
    direct["scope_2"]["equivalence_reference"] = b.fact("1a", f"{b.DOC}:mol:1a", [b.ev(5)])
    direct["scope_bf3"]["charges"]["mass"] = b.fact("2.82 mg", b.q(2.82, "mg"), [b.ev(5)])
    direct["scope_bf3"]["charges"]["mol_percent"] = b.fact("10 mol%", b.q(10, "mol%"), [b.ev(5)])
    direct["scope_bf3"]["equivalence_reference"] = b.fact("1a", f"{b.DOC}:mol:1a", [b.ev(5)])
    direct["scope_dce"]["charges"]["volume"] = b.fact("2.5 mL", b.q(2.5, "mL"), [b.ev(5)])
    scope3a["solvent_conditions"]["composition"] = b.fact("dry 1,2-dichloroethane", ["dry 1,2-dichloroethane"], [b.ev(5)])
    scope3a["solvent_conditions"]["total_volume"] = b.fact("2.5 mL", b.q(2.5, "mL"), [b.ev(5)])
    scope3a["steps"][0] = b.step(
        "reaction", [b.ev(5)], [p["participant_id"] for p in scope3a["participants"]],
        operation="successively add 1a, 2a, BF3·Et2O, and dry DCE; stir",
        temperature=b.q(60, "°C"), time=b.q(12, "h"), atmosphere="argon")
    b.st.write(b.SUBMITTED, bundle)
    stereo = b.st.read(b.RUN / "raw_model/stereo_checks.json")
    for label, mol in (("1c", replacements[f"{b.DOC}:mol:1c"]), ("2f", replacements[f"{b.DOC}:mol:2f"])):
        stereo[label] = b.st.stereo_inspect(b.st.value(mol["smiles"]))
    b.st.write(b.RUN / "raw_model/stereo_checks.json", stereo)


if __name__ == "__main__":
    main()
