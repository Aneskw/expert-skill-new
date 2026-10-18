#!/usr/bin/env python3
"""v3 molecule/reaction contract, bond classification and evidence-preserving import.

No model API is called. Schema checks and graph calculations are not source review.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
VERSION='3.1'
STANDARD='bond-formation-20261016'
def read(p): return json.loads(Path(p).read_text())
def write(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def walk(x):
    if isinstance(x,dict):
        yield x
        for v in x.values(): yield from walk(v)
    elif isinstance(x,list):
        for v in x: yield from walk(v)
def facts(x): return (o for o in walk(x) if 'normalized_value' in o and 'acquisition_method' in o)
def value(f): return f['normalized_value']
def fact(raw=None,normalized=None,evidence=(),method='original_extraction',status=None,verification=None,derivation=None,alternatives=()):
    present=raw is not None or normalized is not None
    return {'raw_value':raw,'normalized_value':normalized,'acquisition_method':method,
     'evidence_ids':list(evidence),'data_status':status or ('pending_review' if present else 'unreported'),
     'verification':verification or ('inferred' if method in ('inference','structure_derivation') and present else 'unreviewed' if present else 'missing'),
     'derivation':derivation,'alternatives':list(alternatives)}
def empty(documents):
    return {'schema_version':'3.1','tool_version':VERSION,'documents':copy.deepcopy(documents),'molecules':[],'reactions':[],'source_context':[]}

def _skeleton(definition,evidence=()):
    """Fill schema boilerplate with explicit unknowns, never chemical defaults."""
    schema=read(ROOT/'references/literature.schema.json')
    def build(node):
        if '$ref' in node:
            name=node['$ref'].split('/')[-1]
            if name=='fact':return fact(evidence=evidence)
            return build(schema['$defs'][name])
        if node.get('type')=='object':return {k:build(v) for k,v in node['properties'].items() if k in node.get('required',[])}
        if node.get('type')=='array':return []
        if 'const' in node:return node['const']
        if 'anyOf' in node or isinstance(node.get('type'),list):return None
        # Metadata enums and IDs require the caller's explicit assignment.
        return None
    return build(schema['$defs'][definition])

def new_molecule(molecule_id,document_id,evidence_ids):
    m=_skeleton('molecule',evidence_ids)
    m.update(molecule_id=molecule_id,identity_status='unresolved',representation_type='unresolved')
    for k in ('chirality_status','geometry_status','representation_status'):
        m['stereochemistry'][k]=fact(normalized='unresolved',evidence=evidence_ids,method='inference',derivation='Scaffold only; source stereochemistry has not been inspected.')
    m['occurrences']=[{'document_id':document_id,'compound_label':fact(evidence=evidence_ids),'evidence_ids':list(evidence_ids)}]
    return m

def new_reaction(reaction_id,document_id,evidence_ids):
    r=_skeleton('reaction',evidence_ids)
    r.update(reaction_id=reaction_id,document_id=document_id,record_type='unresolved',evidence_ids=list(evidence_ids),completeness='partial')
    r['classification']['coverage']=fact(normalized='unresolved',evidence=evidence_ids,method='inference',derivation='Classification has not been performed; no bond-formation labels asserted.')
    r['review']={'status':'needs_review','method':'model','reviewer':None,'checked_fields':[],'unresolved':['Scaffold only: source extraction has not been completed.']}
    return r

def new_participant(participant_id,molecule_id,role,evidence_ids,stage='reaction'):
    p=_skeleton('participant',evidence_ids)
    p.update(participant_id=participant_id,molecule_id=molecule_id,role=role,stage=stage)
    return p

def new_step(step_id,evidence_ids,stage='reaction'):
    s=_skeleton('step',evidence_ids);s.update(step_id=step_id,stage=stage)
    return s

def new_outcome(outcome_id,kind,target_molecule_ids,evidence_ids):
    o=_skeleton('outcome',evidence_ids)
    o.update(outcome_id=outcome_id,kind=kind,target_molecule_ids=list(target_molecule_ids))
    return o

def source_mentions(dataset,terms,context=260):
    """Find literal source-text mentions for cross-region review, never merge facts."""
    import re
    if not terms or any(not t.strip() for t in terms):raise ValueError('Supply nonempty molecule labels/names or source phrases')
    if not 0<=context<=2000:raise ValueError('Context must be between 0 and 2000 characters')
    hits=[]
    for doc in dataset['documents']:
        pages=[e for e in doc['evidence'] if e.get('kind')=='page_text']
        for e in pages:
            text=e.get('text') or ''
            for term in terms:
                pattern=re.escape(term.strip())
                if term.strip().isalnum():pattern=r'(?<!\w)'+pattern+r'(?!\w)'
                for match in re.finditer(pattern,text,re.IGNORECASE):
                    hits.append({'term':term,'document_id':doc['document_id'],'page':e['page'],
                     'evidence_id':e['evidence_id'],'character_range':[match.start(),match.end()],
                     'snippet':text[max(0,match.start()-context):match.end()+context]})
    return {'matches':hits,'limitations':'Literal extracted-text hints only. Inspect original images and associate experiment conditions; identical compound labels do not prove identical experiments. Missing hits do not prove absence from figures, OCR, or the source.'}

def bond_label(a,b,order):
    """Return the exact taxonomy code; atom order follows the standard, not alphabetic order."""
    from rdkit import Chem
    if order not in (1,2,3): raise ValueError('Aromatic/unspecified bonds require source review; no Kekule-based classification')
    table=read(ROOT/'references/bond_taxonomy.json')['categories']
    prefix={1:'',2:'DB-',3:'TB-'}[order]
    for row in table:
        code=row['code']
        if code in ('C-M','Others') or (prefix and not code.startswith(prefix)) or (not prefix and code.startswith(('DB-','TB-'))):continue
        parts=code.removeprefix(prefix).split('-')
        if sorted(parts)==sorted([a,b]): return {'category_number':row['number'],'code':code}
    # Metals are not inferred from the complement of the organic subset: B, Si,
    # halogens, noble gases, etc. must not be mislabeled carbon-metal.
    metal_numbers={3,4,11,12,13,19,20,31,37,38,49,50,55,56,81,82,83,87,88,113,114,115,116}|set(range(21,31))|set(range(39,49))|set(range(57,81))|set(range(89,113))
    other=b if a=='C' else a if b=='C' else None
    if order==1 and other and Chem.GetPeriodicTable().GetAtomicNumber(other) in metal_numbers:
        return {'category_number':23,'code':'C-'+other}
    # The source explicitly uses other-A-B for unlisted formed bonds. Bond order
    # remains in the change, so this fallback never loses double/triple status.
    pair=sorted([a,b],key=lambda s:(s=='H',s))
    return {'category_number':41,'code':'other-'+'-'.join(pair)}

def classify(mapped):
    from rdkit import Chem
    parts=mapped.split('>')
    if len(parts)!=3 or not parts[0] or not parts[2]:raise ValueError('Use atom-mapped reactants>>products; no automatic atom mapping is performed')
    def graph(s):
        mol=Chem.MolFromSmiles(s)
        if mol is None:raise ValueError('Invalid mapped SMILES')
        atoms={};bonds={}
        for a in mol.GetAtoms():
            n=a.GetAtomMapNum()
            if not n or n in atoms:raise ValueError('Every supplied atom needs a unique positive map number on each side')
            if a.GetAtomicNum()==0:raise ValueError('Wildcard/R-group structures cannot be classified as complete molecules')
            atoms[n]={'element':a.GetSymbol(),'hydrogens':a.GetTotalNumHs(),'isotope':a.GetIsotope(),'charge':a.GetFormalCharge()}
        for b in mol.GetBonds():
            bonds[tuple(sorted((b.GetBeginAtom().GetAtomMapNum(),b.GetEndAtom().GetAtomMapNum())))]=b.GetBondTypeAsDouble()
        return atoms,bonds
    ra,rb=graph(parts[0]);pa,pb=graph(parts[2])
    if set(pa)-set(ra):raise ValueError('Product contains atom maps absent from supplied reactants; include the atom source')
    for n in set(ra)&set(pa):
        if (ra[n]['element'],ra[n]['isotope'])!=(pa[n]['element'],pa[n]['isotope']):raise ValueError('Atom map changes element/isotope')
    changes=[];labels={}
    def add(atoms,before,after,kind,hydrogen_count=None):
        c={'change_id':'b'+str(len(changes)+1),'atoms':atoms,'before_order':before,'after_order':after,'change_type':kind}
        if hydrogen_count is not None:c['hydrogen_count']=hydrogen_count
        changes.append(c)
        if kind in ('formed','order_increased','hydrogen_gain'):
            label=bond_label(atoms[0]['element'],atoms[1]['element'],after)
            item=labels.setdefault(label['code'],dict(label,supporting_change_ids=[]))
            item['supporting_change_ids'].append(c['change_id'])
    for pair in sorted(set(rb)|set(pb)):
        before=rb.get(pair,0);after=pb.get(pair,0)
        if before==after:continue
        if before==1.5 or after==1.5:raise ValueError('Changed aromatic bonds require reviewed aromaticity-aware mapping; automatic Kekule classification is disabled')
        kind='formed' if before==0 else 'broken' if after==0 else 'order_increased' if after>before else 'order_decreased'
        add([{'atom_map':n,'element':(pa.get(n) or ra[n])['element']} for n in pair],before,after,kind)
    for n in sorted(set(ra)&set(pa)):
        delta=pa[n]['hydrogens']-ra[n]['hydrogens']
        if delta:
            add([{'atom_map':n,'element':pa[n]['element']},{'atom_map':None,'element':'H'}],0 if delta>0 else 1,1 if delta>0 else 0,'hydrogen_gain' if delta>0 else 'hydrogen_loss',abs(delta))
    return {'standard_id':STANDARD,'bond_changes':changes,'labels':list(labels.values()),
     'coverage':'partial' if any(c['change_type'].startswith('hydrogen_') for c in changes) or set(ra)!=set(pa) else 'mapped_graph_only',
     'limitations':['Calculated from supplied mapping, not verified against the paper.','Implicit-H changes are net gains/losses; proton donor/acceptor and elementary mechanism are not identified.','Unlisted coproducts and spectator changes are outside the supplied mapping.']}

def stereo_inspect(smiles):
    """Graph diagnostics only: unassigned does not mean achiral or racemic.

    Indices are zero-based in the input SMILES. Axial assignment is not certified.
    """
    from rdkit import Chem
    mol=Chem.MolFromSmiles(smiles)
    if mol is None:raise ValueError('Invalid SMILES')
    Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
    potential=[]
    for info in Chem.FindPotentialStereo(mol):
        potential.append({'type':str(info.type),'centered_on':int(info.centeredOn),
                          'specified':str(info.specified),'descriptor':str(info.descriptor)})
    iso=Chem.MolToSmiles(mol,isomericSmiles=True)
    # A graph key is not a sample identity: racemate/unknown/mixture need distinct registry IDs.
    return {'canonical_isomeric_smiles':iso,
            'stereo_graph_key':hashlib.sha256(iso.encode()).hexdigest(),
            'connectivity_smiles':Chem.MolToSmiles(mol,isomericSmiles=False),
            'atom_centers':[{'atom_index':i,'CIP':label} for i,label in Chem.FindMolChiralCenters(mol,includeUnassigned=True)],
            'double_bonds':[{'bond_index':b.GetIdx(),'stereo':str(b.GetStereo())} for b in mol.GetBonds() if b.GetBondType()==Chem.BondType.DOUBLE],
            'potential_stereo':potential,'source_verified':False,
            'limitations':['Graph diagnostics cannot infer unreported configuration, ee, dr or racemic composition.',
                           'Axial/allene, relative and ring cis/trans assignments need source review; toolkit support is limited.']}


def stereo_errors(m):
    errors=[];st=m['stereochemistry'];s=value(m['smiles']);where=m['molecule_id']
    statuses={
      'chirality_status':{'specified','partially_specified','unspecified','racemic','mixture','achiral','unresolved','not_applicable'},
      'geometry_status':{'specified','partially_specified','unspecified','mixture','not_applicable','unresolved'},
      'representation_status':{'encoded','partial','not_encoded','not_applicable','unresolved'}}
    for k,allowed in statuses.items():
        if value(st[k]) is None:errors.append(where+': explicit '+k+' required; use unresolved when uncertain')
        if value(st[k]) is not None and value(st[k]) not in allowed:errors.append(where+': invalid '+k)
    if value(st['chirality_status']) in {'specified','partially_specified'} and not any(value(st[k]) is not None for k in ('atom_centers','axial','assignment')):
        errors.append(where+': specified chirality needs source-linked center/axis/assignment details')
    if value(st['geometry_status']) in {'specified','partially_specified'} and not any(value(st[k]) is not None for k in ('double_bonds','assignment')):
        errors.append(where+': specified geometry needs source-linked E/Z or cis/trans assignment')
    iso=value(st['isomeric_smiles']);key=value(st['identity_key']);diag=None
    if s:
        try:diag=stereo_inspect(s)
        except ValueError:pass  # Main validation reports invalid SMILES.
    if value(st['representation_status'])=='encoded' and not iso:
        errors.append(where+': encoded stereochemistry requires isomeric_smiles')
    if iso:
        try:
            claimed=stereo_inspect(iso)
            if not diag or claimed['canonical_isomeric_smiles']!=diag['canonical_isomeric_smiles']:
                errors.append(where+': isomeric_smiles disagrees with registry SMILES (including stereo)')
        except ValueError:errors.append(where+': invalid isomeric_smiles')
    if diag:
        if key is not None and key!=diag['stereo_graph_key']:errors.append(where+': identity_key must preserve isomeric graph')
        if value(st['chirality_status'])=='achiral' and value(st['assignment']) is None and any(p['type'].startswith('Atom_') and p['specified']=='Unspecified' for p in diag['potential_stereo']):
            errors.append(where+': potential stereocenter cannot be labelled achiral by missing stereo tags')
        if value(st['representation_status'])=='encoded' and any(p['specified']=='Unspecified' for p in diag['potential_stereo']):
            errors.append(where+': encoded claim contains unspecified potential stereochemistry; use partial')
    return errors


def validate(bundle,base=None):
    from jsonschema import Draft202012Validator
    from record_quality import quantity_errors
    errors=[f"schema {'/'.join(map(str,e.path))}: {e.message}" for e in Draft202012Validator(read(ROOT/'references/literature.schema.json')).iter_errors(bundle)]
    if errors:return errors
    docs={d.get('document_id'):d for d in bundle['documents']};evidence={}
    if len(docs)!=len(bundle['documents']):errors.append('duplicate document_id')
    for d in docs.values():
        for e in d.get('evidence',[]):
            eid=e.get('evidence_id')
            if eid in evidence:errors.append('duplicate evidence_id '+str(eid))
            evidence[eid]=e
            if e.get('document_id')!=d.get('document_id') or not isinstance(e.get('page'),int) or not 1<=e['page']<=d.get('pages',0):errors.append('invalid evidence document/page '+str(eid))
        if base:
            paths=[(d.get('source_pdf'),d.get('sha256'))]+[(e['asset'],None) for e in d.get('evidence',[]) if e.get('asset')]
            for relative,digest in paths:
                p=(Path(base)/(relative or '')).resolve()
                if not relative or not p.is_relative_to(Path(base).resolve()) or not p.is_file():errors.append('missing/outside source asset '+str(relative))
                elif digest and sha(p)!=digest:errors.append('source hash mismatch '+str(relative))
    def check_fact(f,allowed,where):
        ev=f['evidence_ids'];present=f['raw_value'] is not None or f['normalized_value'] is not None
        if any(e not in evidence or evidence[e].get('document_id') not in allowed for e in ev):errors.append(where+': unknown or unrelated evidence')
        if present and not ev:errors.append(where+': non-null fact needs evidence')
        if f['data_status'] in ('unreported','not_applicable') and present:errors.append(where+': missing/not-applicable field must be null')
        if f['data_status']=='unreadable' and f['normalized_value'] is not None:errors.append(where+': unreadable normalized value must be null')
        if f['data_status']=='conflict' and (f['normalized_value'] is not None or len(f['alternatives'])<2):errors.append(where+': conflict requires null chosen value and at least two alternatives')
        if not present and f['data_status'] in ('confirmed','pending_review'):errors.append(where+': empty fact needs explicit missing state')
        if present and f['verification']=='missing':errors.append(where+': non-null fact cannot be missing')
        if f['acquisition_method'] in ('inference','structure_derivation') and present:
            if not f['derivation'] or f['verification']!='inferred':errors.append(where+': derived value needs method/derivation/inferred status')
        if f['acquisition_method'] in ('original_extraction','inherited_general_condition','external_database') and f['normalized_value'] is not None and f['raw_value'] is None:errors.append(where+': extraction/inheritance needs raw value')
        for err in quantity_errors(f['normalized_value']):errors.append(where+': '+err)
    mids={m['molecule_id']:m for m in bundle['molecules']}
    if len(mids)!=len(bundle['molecules']):errors.append('duplicate molecule_id')
    for m in bundle['molecules']:
        where=m['molecule_id'];allowed={o['document_id'] for o in m['occurrences']}
        if not allowed or not allowed<=docs.keys():errors.append(where+': molecule requires source occurrences')
        for o in m['occurrences']:
            if not o['evidence_ids'] or any(e not in evidence or evidence[e].get('document_id')!=o['document_id'] for e in o['evidence_ids']):errors.append(where+': invalid occurrence evidence')
        for f in facts(m):check_fact(f,allowed,where)
        errors.extend(stereo_errors(m))
        cas=value(m['cas_numbers'])
        if cas is not None:
            import re
            if not isinstance(cas,list):errors.append(where+': CAS numbers must be a list')
            else:
                for number in cas:
                    if not isinstance(number,str) or not re.fullmatch(r'\d{2,7}-\d{2}-\d',number):errors.append(where+': invalid CAS format')
                    else:
                        digits=number.replace('-','')
                        if sum(int(n)*i for i,n in enumerate(reversed(digits[:-1]),1))%10!=int(digits[-1]):errors.append(where+': CAS checksum mismatch')
        s=value(m['smiles'])
        if m['representation_type']=='smiles':
            if not isinstance(s,str) or not s:errors.append(where+': smiles representation needs SMILES')
            else:
                try:
                    from rdkit import Chem
                    mol=Chem.MolFromSmiles(s)
                    if mol is None or any(a.GetAtomicNum()==0 for a in mol.GetAtoms()):errors.append(where+': invalid or generic SMILES; use template')
                except ImportError:errors.append(where+': RDKit is required to check SMILES')
            if m['generic_template'] is not None:errors.append(where+': choose SMILES or template, not both')
        elif m['representation_type']=='generic_template':
            if m['generic_template'] is None or s is not None:errors.append(where+': template representation requires template and null concrete SMILES')
        elif s is not None or m['generic_template'] is not None:errors.append(where+': unresolved representation cannot claim a complete structure')
        if m['generic_template']:
            for entry in m['generic_template']['instance_mapping']:
                target=mids.get(entry['molecule_id'])
                if not target or target['representation_type']!='smiles':errors.append(where+': template instance must reference a concrete molecule')
        if m['external_properties'] and (m['identity_status']!='structure_checked' or m['representation_type']!='smiles'):
            errors.append(where+': external properties require structure-checked concrete identity')
    reactions={r['reaction_id']:r for r in bundle['reactions']}
    if len(reactions)!=len(bundle['reactions']):errors.append('duplicate reaction_id')
    for r in bundle['reactions']:
        where=r['reaction_id'];did=r['document_id'];allowed={did}|{d['document_id'] for d in docs.values() if d.get('parent_document_id')==did}
        if did not in docs or not where.startswith(did+':'):errors.append(where+': unknown document or non-namespaced reaction_id')
        if not r['evidence_ids'] or any(e not in evidence or evidence[e].get('document_id') not in allowed for e in r['evidence_ids']):errors.append(where+': invalid reaction evidence')
        for f in facts(r):check_fact(f,allowed,where)
        categories=value(r['experiment_categories'])
        if categories is not None and (not isinstance(categories,list) or not set(categories)<={'optimization','substrate_scope','control','scale_up','derivatization','mechanistic_study','other'}):errors.append(where+': invalid experiment categories')
        pids={p['participant_id'] for p in r['participants']};used_mids={p['molecule_id'] for p in r['participants']}
        if len(pids)!=len(r['participants']):errors.append(where+': duplicate participant_id')
        if not used_mids<=mids.keys():errors.append(where+': dangling molecule_id')
        for p in r['participants']:
            if p['molecule_id'] in mids and not any(o['document_id'] in allowed for o in mids[p['molecule_id']]['occurrences']):errors.append(where+': molecule has no occurrence in this paper/SI')
            reference=value(p['equivalence_reference'])
            if reference is not None and reference not in used_mids:errors.append(where+': unknown equivalents/loading denominator')
            if any(value(p['charges'][k]) is not None for k in ('equivalents','mol_percent')) and reference is None:errors.append(where+': equivalents/loading require denominator molecule')
            for k,f in p['charges'].items():
                category={'mol_percent':'percentage'}.get(k,k)
                v=value(f)
                # mol% has its own unit and can exceed 100; do not reuse yield bounds.
                if k=='mol_percent' and v is not None:
                    if not isinstance(v,dict) or v.get('unit')!='mol%':errors.append(where+': loading must retain mol% unit')
                else:
                    for err in quantity_errors(v,category,strict=True):errors.append(where+': '+k+' '+err)
        scale=value(r['reaction_scale'])
        if scale is not None and (not isinstance(scale,dict) or scale.get('reference_molecule_id') not in used_mids):errors.append(where+': reaction scale requires reference molecule')
        for s in r['steps']:
            if not set(s['participant_ids'])<=pids:errors.append(where+': step has unknown participant')
            for k in ('temperature','time'):
                for err in quantity_errors(value(s[k]),k,strict=True):errors.append(where+': '+k+' '+err)
        for o in r['outcomes']:
            if not set(o['target_molecule_ids'])<=used_mids:errors.append(where+': outcome has dangling target')
            if o['kind'] in ('yield','product_amount','selectivity','recovery','conversion') and not o['target_molecule_ids']:errors.append(where+': quantitative outcome requires target')
            if o['kind']=='yield' and any(not any(p['molecule_id']==t and p['role']=='product' for p in r['participants']) for t in o['target_molecule_ids']):errors.append(where+': yield target must be a product')
            if o['kind']=='recovery' and any(not any(p['molecule_id']==t and p['role']!='product' for p in r['participants']) for t in o['target_molecule_ids']):errors.append(where+': recovery target must be input material')
            if o['kind'] in ('yield','conversion','recovery'):
                for err in quantity_errors(value(o['value']),o['kind'],strict=True):errors.append(where+': '+o['kind']+' '+err)
        for sid in r['general_scheme_ids']:
            if sid not in reactions or reactions[sid]['record_type']!='generic_scheme':errors.append(where+': general_scheme_id must reference a generic scheme')
        c=r['classification'];labels=value(c['labels']);coverage=value(c['coverage'])
        if coverage not in ('complete','partial','unresolved','not_applicable'):errors.append(where+': invalid classification coverage')
        changes={value(f).get('change_id'):value(f) for f in c['bond_changes'] if isinstance(value(f),dict)}
        if len(changes)!=len(c['bond_changes']):errors.append(where+': duplicate/malformed bond change')
        if labels is not None and not isinstance(labels,list):errors.append(where+': classification labels must be a list')
        if labels and coverage in ('unresolved','not_applicable'):errors.append(where+': unresolved/not-applicable classification cannot assert labels')
        if coverage=='not_applicable' and (c['bond_changes'] or labels not in ([],None)):errors.append(where+': no applicable transformation must not invent bond changes')
        expected={}
        try:
            for change_id,ch in changes.items():
                before,after=ch['before_order'],ch['after_order'];kind=ch['change_type']
                if kind not in ('formed','broken','order_increased','order_decreased','hydrogen_gain','hydrogen_loss'):raise ValueError('unknown bond change type')
                if kind in ('formed','hydrogen_gain') and (before!=0 or after not in (1,2,3)):raise ValueError('formed bond requires absent-to-present change')
                if kind=='order_increased' and not 0<before<after:raise ValueError('bond-order increase must increase order')
                if kind=='order_decreased' and not before>after>0:raise ValueError('bond-order decrease must decrease positive order')
                if kind in ('broken','hydrogen_loss') and not before>after==0:raise ValueError('broken bond requires present-to-absent change')
                if kind in ('formed','hydrogen_gain','order_increased'):
                    lab=bond_label(ch['atoms'][0]['element'],ch['atoms'][1]['element'],after)
                    expected.setdefault(lab['code'],dict(lab,supporting_change_ids=[]))['supporting_change_ids'].append(change_id)
            actual={l['code']:l for l in labels or []}
            if labels is not None and (len(actual)!=len(labels) or actual!=expected):errors.append(where+': taxonomy labels must match all and only formed/promoted bonds')
            mapped=value(c['mapped_reaction_smiles'])
            if mapped:
                calculated=classify(mapped)
                if [value(f) for f in c['bond_changes']]!=calculated['bond_changes']:errors.append(where+': bond changes do not match supplied atom mapping')
                from rdkit import Chem
                def components(s):
                    mol=Chem.MolFromSmiles(s)
                    for a in mol.GetAtoms():a.SetAtomMapNum(0)
                    return set(Chem.MolToSmiles(f,isomericSmiles=True) for f in Chem.GetMolFrags(mol,asMols=True))
                for side,roles in [(0,{'reactant','reagent','solvent','additive'}),(2,{'product'})]:
                    known=set()
                    for p in r['participants']:
                        m=mids.get(p['molecule_id'],{})
                        if p['role'] in roles and m.get('representation_type')=='smiles':known|=components(value(m['smiles']))
                    if not components(mapped.split('>')[side])<=known:errors.append(where+': mapped structures do not match registered reaction participants')
                # Implicit H mapping and omitted partners cannot be called complete.
                if coverage=='complete' and calculated['coverage']!='mapped_graph_only':errors.append(where+': incomplete hydrogen/partner mapping cannot be complete')
        except (ValueError,KeyError,TypeError,ImportError,RuntimeError) as exc:errors.append(where+': classification '+str(exc))
        if coverage=='complete' and not value(c['mapped_reaction_smiles']):errors.append(where+': complete classification requires atom mapping')
        if r['review']['status']=='source_reviewed' and r['review']['method'] not in ('human','assistant_visual'):errors.append(where+': source review is not model consensus')
        if r['completeness']=='complete' and (r['review']['unresolved'] or not r['steps'] or r['record_type']!='specific_experiment'):errors.append(where+': incomplete experiment cannot claim complete')
    return sorted(set(errors))

def quarantine(item,is_molecule=False):
    item=copy.deepcopy(item)
    for f in facts(item):
        if f['data_status']=='confirmed':f['data_status']='pending_review'
        if f['verification']=='source_checked':f['verification']='unreviewed'
    if is_molecule:
        if item['identity_status']=='structure_checked':item['identity_status']='source_identified'
    else:item['review']={'status':'needs_review','method':'model','reviewer':None,'checked_fields':[],'unresolved':item['review']['unresolved']}
    return item

def legacy_fact(f):
    v=value(f)
    return {'value':copy.deepcopy(v),'evidence_ids':list(f['evidence_ids']),
     'verification':'missing' if v is None else f['verification']}

def to_legacy(r,molecules):
    """Compatibility view only. v3 remains the full source of raw values and identities."""
    from reaction_tools import smiles_ok
    kind={'generic_scheme':'general_scheme','proposed_mechanism':'mechanism','literature_example':'literature_example','unresolved':'unresolved'}.get(r['record_type'],'reported_experiment')
    categories=value(r['experiment_categories']) or []
    if kind=='reported_experiment':
        if 'optimization' in categories:kind='optimization_entry'
        elif 'substrate_scope' in categories:kind='substrate_scope_entry'
    result={'record_id':r['reaction_id'],'document_id':r['document_id'],'source_reaction_id':r['reaction_id'].split(':',1)[1],
     'record_kind':kind,'status':'needs_review','completeness':'partial','evidence_ids':r['evidence_ids'],
     'participants':[],'steps':[],'outcomes':[],'notes':[],
     'review':{'method':'model','reviewer':None,'checked_fields':[],'unresolved':r['review']['unresolved']}}
    def field(name,f,targets=()):return {'name':name,'fact':legacy_fact(f),'target_entity_ids':list(targets)}
    for p in r['participants']:
        m=molecules[p['molecule_id']];occ=next((o for o in m['occurrences'] if o['document_id']==r['document_id']),m['occurrences'][0])
        result['participants'].append({'entity_id':p['participant_id'],
         'role':p['role'] if p['role'] in ('reactant','product','catalyst','reagent','solvent') else 'unspecified_agent',
         'label':legacy_fact(occ['compound_label']),'name':legacy_fact(m['names']['original']),
         'smiles':legacy_fact(m['smiles']),'smiles_parseable':smiles_ok(value(m['smiles'])),'structure_verified':False,
         'amounts':[field('mol_percent' if k=='mol_percent' else k,f) for k,f in p['charges'].items()]})
    for s in r['steps']:
        result['steps'].append({'step_id':s['step_id'],'participant_ids':s['participant_ids'],
         'conditions':[field(k,s[k]) for k in ('operation','temperature','time','atmosphere','pressure')]})
    for o in r['outcomes']:
        targets=[p['participant_id'] for p in r['participants'] if p['molecule_id'] in o['target_molecule_ids'] and (p['role']=='product' if o['kind'] in ('yield','product_amount') else True)]
        result['outcomes'] += [field(o['kind'],o['value'],targets),field(o['kind']+'_measurement_basis',o['measurement_basis'],targets)]
    result['notes']=[legacy_fact(f) for f in r['notes']]
    return result

def merge(dataset_path,bundle_path):
    from reaction_tools import validate as legacy_validate,normalize_record,report
    path=Path(dataset_path);base=path.resolve().parent;data=read(path);submitted=read(bundle_path)
    trial=copy.deepcopy(data);old=data.get('standard',empty(data['documents']))
    incoming=copy.deepcopy(submitted)
    # Do not accept submitted evidence additions/replacement through machine merge.
    if incoming.get('documents')!=data['documents']:raise ValueError('Submitted documents/evidence must exactly match the prepared dataset')
    errors=validate(incoming,base)
    if errors:raise ValueError('Standard merge refused: '+json.dumps(errors,ensure_ascii=False))
    molecules={m['molecule_id']:m for m in old['molecules']};reactions={r['reaction_id']:r for r in old['reactions']}
    for m in incoming['molecules']:
        existing=molecules.get(m['molecule_id'])
        if existing==m:continue
        if existing and any(f['verification']=='source_checked' for f in facts(existing)):raise ValueError('Cannot overwrite source-reviewed molecule')
        molecules[m['molecule_id']]=quarantine(m,True)
    for r in incoming['reactions']:
        rid=r['reaction_id']
        if reactions.get(rid,{}).get('review',{}).get('status')=='source_reviewed':raise ValueError('Cannot overwrite source-reviewed reaction')
        if any(e['record_id']==rid and e.get('status')=='reviewed' for e in data['experiments']):raise ValueError('Cannot overwrite legacy source-reviewed reaction')
        reactions[rid]=quarantine(r)
    standard=empty(data['documents']);standard['molecules']=list(molecules.values());standard['reactions']=list(reactions.values())
    standard['source_context']=old.get('source_context',[])+[c for c in incoming['source_context'] if c not in old.get('source_context',[])]
    errors=validate(standard,base)
    if errors:raise ValueError('Standard merged bundle invalid: '+json.dumps(errors,ensure_ascii=False))
    indexed={r['record_id']:r for r in trial['experiments']}
    for r in incoming['reactions']:
        legacy=to_legacy(reactions[r['reaction_id']],molecules);normalize_record(legacy);indexed[legacy['record_id']]=legacy
    trial['experiments']=list(indexed.values());trial['standard']=standard;trial['tool_version']=VERSION
    errors=legacy_validate(trial,base,strict_record_ids={r['reaction_id'] for r in incoming['reactions']})
    if errors:raise ValueError('Compatibility view invalid: '+json.dumps(errors,ensure_ascii=False))
    write(path,trial);write(base/'standard_records.json',standard);render(standard,base);report(base)
    return {'merged_reactions':len(incoming['reactions']),'molecules':len(molecules),'chemical_accuracy_verified':False}

def render(bundle,out):
    """Human-readable complete records plus multi-label taxonomy index; no copy per class."""
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    mols={m['molecule_id']:m for m in bundle['molecules']};ev={e['evidence_id']:e for d in bundle['documents'] for e in d.get('evidence',[])}
    index={};pages=['# 按规范提取的反应记录','','数据来自文献；来源核对不等于实验复现。每条记录中的缺失字段、推断与冲突均保留。','']
    def display(f):
        return json.dumps(f,ensure_ascii=False)
    for r in bundle['reactions']:
        labels=value(r['classification']['labels'])
        for l in labels or []:index.setdefault(l['code'],[]).append(r['reaction_id'])
        if not labels:index.setdefault('未分类或不适用',[]).append(r['reaction_id'])
        lines=['# '+r['reaction_id'],'',f"记录类型：{r['record_type']}；复核状态：{r['review']['status']}；完整性：{r['completeness']}",'',
         '分类：'+', '.join(l['code'] for l in labels or [])+'；覆盖：'+str(value(r['classification']['coverage'])),'']
        def emit(x,path=''):
            if isinstance(x,dict) and 'normalized_value' in x and 'acquisition_method' in x:
                lines.append('- `'+path+'`：'+display(x))
            elif isinstance(x,dict):
                for k,v in x.items():emit(v,path+'.'+k if path else k)
            elif isinstance(x,list):
                for i,v in enumerate(x):emit(v,path+'['+str(i)+']')
            elif x is not None:lines.append('- `'+path+'`：'+str(x))
        emit(r)
        for mid in dict.fromkeys(p['molecule_id'] for p in r['participants']):
            lines+=['','## 分子 '+mid,''];emit(mols[mid],mid)
        cited={eid for f in facts(r) for eid in f['evidence_ids']}|set(r['evidence_ids'])
        lines+=['','## 原文定位','']
        for eid in sorted(cited):
            if eid in ev:
                e=ev[eid];lines.append(f"- {eid}：PDF 第 {e['page']} 页"+(f"；[原图]({e['asset']})" if e.get('asset') else ''))
        stem=r['reaction_id'].replace(':','__')
        write(out/'records'/(stem+'.json'),r)
        # Link assets relative to records/ in per-record Markdown.
        (out/'records'/(stem+'.md')).write_text('\n'.join(lines).replace('](assets/','](../assets/')+'\n')
        pages+=['## '+r['reaction_id'],'','[完整记录](records/'+stem+'.md)','']+lines[2:6]+['']
    write(out/'molecules.json',bundle['molecules']);write(out/'classification_index.json',index)
    pages+=['## 分类索引','','同一反应可进入多个形成键类别，原记录只保留一份。','']
    for label,ids in index.items():pages.append('- '+label+'：'+', '.join(ids))
    (out/'规范记录预览.md').write_text('\n'.join(pages)+'\n')
    return index

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('init');a.add_argument('--dataset',required=True);a.add_argument('--out',required=True)
    a=sub.add_parser('validate');a.add_argument('--bundle',required=True)
    a=sub.add_parser('merge');a.add_argument('--dataset',required=True);a.add_argument('--bundle',required=True)
    a=sub.add_parser('classify');a.add_argument('--mapped-reaction',required=True)
    a=sub.add_parser('mentions');a.add_argument('--dataset',required=True);a.add_argument('--term',action='append',required=True);a.add_argument('--context',type=int,default=260)
    a=sub.add_parser('render');a.add_argument('--bundle',required=True);a.add_argument('--out',required=True)
    a=sub.add_parser('stereo');a.add_argument('--smiles',required=True)
    args=p.parse_args(argv)
    try:
        if args.command=='init':
            if Path(args.out).exists():raise ValueError('Refuse to overwrite existing bundle')
            write(args.out,empty(read(args.dataset)['documents']))
        elif args.command=='validate':
            errors=validate(read(args.bundle),Path(args.bundle).resolve().parent)
            print(json.dumps({'valid':not errors,'errors':errors,'chemical_accuracy_verified':False},ensure_ascii=False,indent=2))
            if errors:return 1
        elif args.command=='merge':print(json.dumps(merge(args.dataset,args.bundle),ensure_ascii=False))
        elif args.command=='classify':print(json.dumps(classify(args.mapped_reaction),ensure_ascii=False,indent=2))
        elif args.command=='stereo':print(json.dumps(stereo_inspect(args.smiles),ensure_ascii=False,indent=2))
        elif args.command=='mentions':print(json.dumps(source_mentions(read(args.dataset),args.term,args.context),ensure_ascii=False,indent=2))
        elif args.command=='render':
            b=read(args.bundle);errors=validate(b,Path(args.bundle).resolve().parent)
            if errors:raise ValueError('Refuse to render invalid bundle: '+str(errors))
            if Path(args.out).resolve()!=Path(args.bundle).resolve().parent:raise ValueError('Render beside bundle to preserve relative source links')
            print(json.dumps(render(b,args.out),ensure_ascii=False))
    except (ValueError,KeyError,ImportError,OSError) as e:print('ERROR: '+str(e),file=sys.stderr);return 1
    return 0
if __name__=='__main__':sys.exit(main())
