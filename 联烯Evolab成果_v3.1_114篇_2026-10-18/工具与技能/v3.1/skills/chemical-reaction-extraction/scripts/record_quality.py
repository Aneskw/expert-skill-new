"""Deterministic record checks; these do not certify source support."""
import math
import re

VERSION = '3.1'
UNIT_ALIASES = {'days': 'd', 'day': 'd', 'hours': 'h', 'hour': 'h',
                'minutes': 'min', 'seconds': 's', 'degC': '°C', 'Celsius': '°C',
                '℃': '°C', 'eq': 'equiv', 'eq.': 'equiv', 'equivalents': 'equiv',
                'μmol': 'umol', 'µmol': 'umol', 'μL': 'uL', 'µL': 'uL'}
NONNEGATIVE = {'%', 'mol%', 'equiv', 'd', 'h', 'min', 's', 'mg', 'g', 'kg',
               'ug', 'mmol', 'mol', 'umol', 'nmol', 'mL', 'L', 'uL', 'M', 'mM', 'K'}
UNITS = NONNEGATIVE | {'°C', '°F'}
NUMBER = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)'
SCALAR = re.compile(rf'^\s*(?P<qualifier><=|>=|≤|≥|<|>|~|≈)?\s*(?P<number>{NUMBER})\s*(?P<unit>[^\d\s]+)\s*$')
RANGE = re.compile(rf'^\s*(?P<minimum>{NUMBER})\s*(?:–|—|to|-)\s*(?P<maximum>{NUMBER})\s*(?P<unit>[^\d\s]+)\s*$')
QUALIFIERS = {'<': 'lt', '>': 'gt', '<=': 'le', '≤': 'le', '>=': 'ge', '≥': 'ge', '~': 'approx', '≈': 'approx'}


def category(name):
    name = name.casefold().replace('-', '_')
    if 'measurement' in name or name.endswith('_basis'):
        return 'measurement_basis'
    if 'recover' in name or name.startswith(('starting_material_yield', 'yield_of_starting_material', 'unreacted_starting_material_yield')):
        return 'recovery'
    if name == 'conversion' or name.endswith('_conversion'):
        return 'conversion'
    if 'yield' in name:
        return 'yield'
    if 'time' in name.split('_'):
        return 'time'
    if 'temperature' in name or name == 'warm_to':
        return 'temperature'
    if 'mass' in name:
        return 'mass'
    if 'volume' in name:
        return 'volume'
    if 'equiv' in name:
        return 'equivalents'
    if name in ('amount', 'reaction_scale'):
        return 'amount'
    if 'concentration' in name:
        return 'concentration'
    if name in ('ee', 'de', 'conversion', 'enantiomeric_excess') or name.endswith(('_ee', '_de')):
        return 'percentage'
    if name == 'operation' or name.endswith('_operation'):
        return 'operation'
    if name == 'solvent' or name.endswith('_solvent'):
        return 'solvent'
    if {'ratio', 'ratios'}.intersection(name.split('_')) or 'selectivity' in name:
        return 'selectivity'
    return name if name in ('atmosphere', 'workup', 'purification', 'operation', 'solvent') else 'other'


def normalize(value):
    """Parse explicit units only; preserve the original expression."""
    if isinstance(value, dict):
        result = dict(value)
        if 'unit' in result:
            normalized = UNIT_ALIASES.get(result['unit'], result['unit'])
            if normalized != result['unit']:
                result.setdefault('original_unit', result['unit'])
            result['unit'] = normalized
        return result
    if not isinstance(value, str):
        return value
    match = RANGE.fullmatch(value)
    if match and UNIT_ALIASES.get(match['unit'], match['unit']) in UNITS:
        return {'minimum': float(match['minimum']), 'maximum': float(match['maximum']),
                'unit': UNIT_ALIASES.get(match['unit'], match['unit']), 'raw_text': value}
    match = SCALAR.fullmatch(value)
    if match and UNIT_ALIASES.get(match['unit'], match['unit']) in UNITS:
        result = {'number': float(match['number']), 'unit': UNIT_ALIASES.get(match['unit'], match['unit']), 'raw_text': value}
        if match['qualifier']:
            result['qualifier'] = QUALIFIERS[match['qualifier']]
        return result
    return value


def quantity_errors(value, field_category=None, strict=False):
    original = value
    value = normalize(value)
    errors = []
    if isinstance(value, dict) and 'unit' in value:
        unit = value['unit']
        scalar = 'number' in value
        ranged = 'minimum' in value or 'maximum' in value
        if scalar == ranged:
            return ['quantity needs either number or minimum/maximum']
        numbers = [value.get('number')] if scalar else [value.get('minimum'), value.get('maximum')]
        if any(isinstance(n, bool) or not isinstance(n, (int, float)) or not math.isfinite(n) for n in numbers):
            return ['quantity requires finite numeric values']
        if ranged and numbers[0] > numbers[1]:
            errors.append('quantity range reversed')
        if unit in NONNEGATIVE and any(n < 0 for n in numbers):
            errors.append('negative non-temperature quantity')
        # Excess catalyst loading can legitimately exceed 100 mol%; yields cannot.
        if unit == '%' and any(not 0 <= n <= 100 for n in numbers):
            errors.append('percentage outside 0..100')
        if value.get('qualifier') not in (None, 'lt', 'le', 'gt', 'ge', 'approx'):
            errors.append('unknown quantity qualifier')
        expected = {'time': {'d', 'h', 'min', 's'}, 'temperature': {'°C', '°F', 'K'},
                    'yield': {'%'}, 'recovery': {'%'}, 'conversion': {'%'}, 'percentage': {'%'}, 'mass': {'ug', 'mg', 'g', 'kg'},
                    'volume': {'uL', 'mL', 'L'}, 'equivalents': {'equiv'}}
        if field_category in expected and unit not in expected[field_category]:
            errors.append('unit does not match field category')
        if strict and unit not in UNITS:
            errors.append('unsupported unit; preserve as unresolved text until supported')
    elif strict and original is not None and field_category in ('time', 'temperature', 'yield', 'recovery', 'conversion', 'percentage', 'mass', 'volume', 'equivalents'):
        # Qualitative statements are meaningful; never convert them to zero.
        if not isinstance(original, str) or re.search(r'\d', original):
            errors.append('numeric field needs explicit structured quantity and unit')
    return errors


def fields(record):
    for participant in record.get('participants', []):
        yield from participant.get('amounts', [])
    for step in record.get('steps', []):
        yield from step.get('conditions', [])
    yield from record.get('outcomes', [])


def normalize_record(record):
    for field in fields(record):
        field['category'] = category(field['name'])
        if field['category'] in ('time', 'temperature', 'yield', 'recovery', 'conversion', 'percentage', 'mass', 'volume', 'equivalents', 'amount', 'concentration'):
            field['fact']['value'] = normalize(field['fact']['value'])
    return record


def record_errors(record, strict=False):
    errors = []
    ids = {p['entity_id'] for p in record.get('participants', [])}
    for field in fields(record):
        inferred = category(field['name'])
        actual = field.get('category', inferred)
        if strict and inferred != 'other' and actual != inferred:
            errors.append('field category disagrees with name')
        if any(identity not in ids for identity in field.get('target_entity_ids', [])):
            errors.append('field refers to unknown target entity')
        errors.extend(quantity_errors(field['fact']['value'], actual, strict))
    steps = [s['step_id'] for s in record.get('steps', [])]
    if len(steps) != len(set(steps)):
        errors.append('duplicate step_id')
    if record.get('completeness') == 'complete':
        if record.get('record_kind') not in ('reported_experiment', 'optimization_entry', 'substrate_scope_entry'):
            errors.append('complete requires an experimental record kind')
        roles = {p['role'] for p in record.get('participants', [])}
        if not {'reactant', 'product'} <= roles or not steps:
            errors.append('complete requires reactant, product and experimental steps')
        if record.get('review', {}).get('unresolved'):
            errors.append('complete record still has unresolved items')
        all_fields = list(fields(record))
        supplied = {category(f['name']) for f in all_fields if f['fact']['value'] is not None}
        if not {'temperature', 'time', 'yield'} <= supplied:
            errors.append('complete requires explicit temperature, time and outcome; otherwise retain partial')
    return errors


def warnings(record):
    result = []
    for step in record.get('steps', []):
        conditions = step.get('conditions', [])
        temperatures = [f['fact']['value'] for f in conditions if category(f['name']) == 'temperature']
        for field in conditions:
            match = re.search(r'time_at_minus_(\d+)(?:_|$)', field['name'])
            if match and not any(isinstance(v, dict) and v.get('number') == -int(match[1]) and v.get('unit') == '°C' for v in temperatures):
                result.append('stage temperature appears only in time field name; add paired temperature evidence: ' + field['name'])
            if field['name'] == 'time_at_room_temperature' and not any(isinstance(v, str) and v.casefold() in ('r.t.', 'rt', 'room temperature') for v in temperatures):
                result.append('room-temperature hold needs its own paired temperature fact')
    for field in fields(record):
        if field.get('category') and category(field['name']) != field['category']:
            result.append('stored field category differs from current classification: ' + field['name'])
    # Advisory only: no chemistry role is inferred from a remark.
    quantities = [str(f['fact']['value']) for p in record.get('participants', []) for f in p.get('amounts', [])]
    for note in record.get('notes', []):
        text = note.get('value')
        if isinstance(text, str) and re.search(r'no\s+(?:yield\s+)?improvement|did\s+not\s+improve', text, re.I):
            mentioned = re.findall(r'(\d+(?:\.\d+)?)\s*(equiv|mmol|mg|mL)\b', text)
            existing = [normalize(f['fact']['value']) for p in record.get('participants', []) for f in p.get('amounts', [])]
            if mentioned and all(any(isinstance(v, dict) and v.get('number') == float(number) and v.get('unit') == unit for v in existing) for number, unit in mentioned):
                # Discussion of already represented loadings is not a missing addition.
                continue
        if isinstance(text, str) and re.search(r'\d\s*(?:equiv|mmol|mol%|mg|mL)\b', text) and not quantities:
            result.append('quantity appears in notes but no participant amounts are structured')
        elif isinstance(text, str) and re.search(r'\d\s*equiv\s+\S+', text):
            result.append('check that note reagent and equivalents also appear in participants/amounts')
    if record.get('record_kind') in ('reported_experiment', 'optimization_entry', 'substrate_scope_entry') and not record.get('participants'):
        result.append('experimental record has no structured participants')
    if record.get('record_kind') in ('reported_experiment', 'optimization_entry', 'substrate_scope_entry'):
        participants = record.get('participants', [])
        reactants = [p for p in participants if p.get('role') == 'reactant']
        products = [p for p in participants if p.get('role') == 'product']
        outcomes = record.get('outcomes', [])
        active = [f for f in outcomes if f.get('fact', {}).get('value') is not None]
        if not reactants:
            result.append('reaction_identity: no reactant; check whether this is characterization or a non-synthetic observation')
        if not record.get('steps'):
            result.append('missing_steps: recover reported procedure or explicitly retain missing conditions')
        for f in active:
            kind = category(f['name'])
            if kind in ('yield', 'recovery', 'conversion', 'selectivity', 'percentage') and not f.get('target_entity_ids'):
                result.append('outcome_target: explicit entity association missing for ' + f['name'])
            if kind in ('yield', 'recovery', 'conversion'):
                bases = [b for b in active if category(b['name']) == 'measurement_basis']
                targets = set(f.get('target_entity_ids', []))
                if not any(not b.get('target_entity_ids') or targets.intersection(b['target_entity_ids']) for b in bases):
                    result.append('measurement_basis: missing or unassigned method for ' + f['name'] + '; do not assume isolated yield')
            value = normalize(f['fact']['value'])
            if kind == 'yield' and not products and isinstance(value, dict) and value.get('number', value.get('minimum', 0)) > 0:
                result.append('product_identity: positive yield without a product participant')
        # Do not demand a fabricated product for failed/no-reaction experiments.
        primary = reactants + products
        if primary and all(p.get('smiles', {}).get('value') is None for p in primary):
            result.append('structure_gap: all reactant/product SMILES missing; inspect drawings and explain ambiguity without inventing structures')
        if not reactants and any(re.search(r'peak_list|molecular_weight|spectral|chemical_shift|relative_abundance', f['name'], re.I) for f in active):
            result.append('characterization_only: spectra/fragment abundance are not synthetic reaction experiments; preserve as source context')
    return sorted(set(result))
