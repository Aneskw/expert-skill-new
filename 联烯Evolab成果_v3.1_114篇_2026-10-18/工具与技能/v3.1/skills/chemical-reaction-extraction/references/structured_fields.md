# Structured quantities and participants

Read this when writing quantities, table footnotes or product-specific outcomes.

Preserve `name` as the detailed field name. `merge` adds a standard `category` (time, temperature, yield, percentage, mass, volume, equivalents, etc.) without collapsing NMR versus isolated yields or combining stages. Use `target_entity_ids` to associate an outcome with its product participant; these IDs must exist in the record. A combined product yield can cite several target IDs. Keep the measurement basis as a separate field with evidence.

Outcome categories are distinct: `yield` for produced material, `recovery` for recovered starting material/reagent, `conversion` for consumption of starting material. The phrase “yield of starting material” in a footnote belongs to recovery when it reports retained starting material. Preserve the author's detailed field name, numeric value, measurement method and entity association; do not compare recovery percentages as product yields. For an intermediate made in one stage then consumed later, retain its stage yield rather than forcing a global product-only role rule.

Use explicit quantities:

- Exact: `{"number": 82, "unit": "%"}`.
- Upper bound: `{"number": 5, "unit": "%", "qualifier": "lt", "raw_text": "<5%"}`.
- Approximate: `{"number": 10, "unit": "%", "qualifier": "approx", "raw_text": "~10%"}`.
- Range: `{"minimum": 65, "maximum": 70, "unit": "°C", "raw_text": "65–70 °C"}`.
- Qualitative: preserve `"trace"`, `"not detected"`, `"r.t."`, `"reflux"`, or equivalent original text. Do not turn them into zero or an assumed temperature.
- Ambiguous numeric statement: keep the original expression in notes, set the relevant fact to null/missing and explain the unresolved unit or assignment. Do not force unsupported units or multi-stage programs into a single number.

Supported units include d, h, min, s, °C, °F, K, %, mol%, equiv, ug, mg, g, kg, umol, mmol, mol, nmol, uL, mL, L, M and mM. The helper normalizes explicit aliases such as days→d, degC→°C, eq→equiv. It parses only explicit units; a table header may supply the unit if its evidence is cited by the model. Qualifiers never remove the original text.

Every reagent or additive explicitly introduced by a footnote belongs in `participants`, with its quantity in `amounts`, plus the footnote evidence. Notes may explain exceptions but must not be the only place a known reagent or its known amount appears. Preserve an uncertain chemical role as `unspecified_agent`; do not infer a catalyst role merely from low loading.

Represent a multistage temperature program as ordered `steps`, each with its own `temperature` and `reaction_time`/`addition_time` facts. Preserve separate addition and hold periods. Include explicitly timed ramps using `temperature_start`, `temperature_end` and `transition_time`; leave unspecified ramp times missing. Room temperature remains qualitative. Keep reaction-solvent volumes, aldehyde-solution volumes and reagent-solution carriers distinguishable; do not assume their sum is an author-reported total volume.

If a scheme reports a loading but the prose reports a different charge, preserve the scheme loading and prose amount independently. Record the conflict in evidence-linked notes and `review.unresolved`. Calculating an equivalent from charge ratios requires an explicit denominator and `verification=inferred`; a crude intermediate is not automatically quantitatively present at its precursor charge. Do not silently repair the author's numbers. Likewise, an overall or two-step isolated yield must retain its starting-material-to-product span and must not be relabeled as the last-step yield.

Examples of checks to perform using `audit`: an added reagent quantity stranded in notes, an experimental record with no participants. Warnings are advisory and need source-based resolution, not blind edits. Log either the correction or why the information remains unresolved.

The v2.3 advisory checks also expose unassigned outcomes, missing measurement methods, reaction-identity problems, absent steps and missing primary structures. They do not prove that a field was reported in the source. In particular, failed experiments may have no product, standalone spectra should remain source context rather than synthetic experiments, and unavailable SI must remain a limitation. Never fill these gaps using chemical convention.

Use `yield_measurement_basis` (or another name ending `_basis`) with the same `target_entity_ids` as its yield when products have different methods. A shared basis can remain untargeted only when the source explicitly applies it to the entire record. A method stated in a field name is not a substitute for the explicit fact. The Markdown includes entity IDs, per-field evidence IDs and verification states so a reviewer can trace each number rather than only the record as a whole.

An experimental `complete` record requires reactant and product participants, steps, explicit time/temperature/outcome, no unresolved items, and the existing SI/not-applicable check. If these criteria do not fit the source, retain `partial`; completeness is not a claim that chemistry has been independently verified.

For structures, separate recognizing a molecule from verifying its drawing. If a full structure is legible, attempt an evidence-linked candidate SMILES using available visual tools or provided chemistry candidates; retain `structure_verified=false`. If the drawing is a generic R-group or ambiguous, leave SMILES null and record a reason. Neither model confidence nor RDKit parsing establishes source correctness.
