#!/usr/bin/env python3
"""Evidence-linked reaction extraction, tools for Evolab; the agent selects the extraction workflow."""
from pathlib import Path
import argparse,copy,hashlib,json,math,re,shutil,sqlite3,sys
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(Path(__file__).resolve().parent))
from record_quality import VERSION, category, normalize_record, quantity_errors, record_errors, warnings

def read(p):return json.loads(Path(p).read_text())
def write(p,d):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False));tmp.replace(p)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fact(value,ev=(),verification='unreviewed'):
 return {'value':value,'evidence_ids':list(ev),'verification':'missing' if value is None else verification}
def walk(x):
 if isinstance(x,dict):
  yield x
  for v in x.values():yield from walk(v)
 elif isinstance(x,list):
  for v in x:yield from walk(v)
def smiles_ok(s):
 if not s:return None
 try:
  from rdkit import Chem,RDLogger
  RDLogger.DisableLog('rdApp.*');return Chem.MolFromSmiles(s) is not None
 except ImportError:return None

def add_asset(src,out,docid,name):
 dest=out/'assets'/docid/name;dest.parent.mkdir(parents=True,exist_ok=True)
 if not dest.exists() or sha(src)!=sha(dest):shutil.copy2(src,dest)
 return dest.relative_to(out).as_posix()

def ingest(pdf,folder,out,role='main_text',parent=None,render=False):
 from pypdf import PdfReader
 pdf=Path(pdf);digest=sha(pdf);docid='doc_'+digest[:16];r=PdfReader(pdf);evidence=[]
 doc={'document_id':docid,'sha256':digest,'source_name':pdf.name,'source_pdf':add_asset(pdf,out,docid,'source.pdf'),'role':role,'parent_document_id':parent,'pages':len(r.pages),'title':None,'doi':None,'supplementary_material':{'status':'unknown','note':'SI availability is not established from filename.'},'evidence':evidence,'coverage':{'pages_processed':len(r.pages),'semantic_extraction_complete':False,'note':'Page/block enumeration is not proof of complete reaction extraction.'}}
 for i,p in enumerate(r.pages):
  text=p.extract_text() or '';evidence.append({'evidence_id':f'{docid}:page{i+1}:text','document_id':docid,'page':i+1,'kind':'page_text','text':text,'asset':None,'bbox':None,'bbox_units':None,'provenance':'pypdf'})
  if i==0:
   m=re.search(r'10\.\d{4,9}/[^\s<>]+',text);doc['doi']=m.group(0).rstrip('.,;') if m else None
 if render:
  import pypdfium2 as pdfium
  engine=pdfium.PdfDocument(str(pdf))
  try:
   for i in range(len(engine)):
    page=engine[i]
    try:
     bitmap=page.render(scale=2)
     try:
      dest=out/'assets'/docid/f'page{i+1}.jpg';bitmap.to_pil().save(dest,quality=90)
     finally:bitmap.close()
    finally:page.close()
    evidence.append({'evidence_id':f'{docid}:page{i+1}:image','document_id':docid,'page':i+1,'kind':'page_image','text':None,'asset':dest.relative_to(out).as_posix(),'bbox':None,'bbox_units':None,'provenance':'rendered_source_pdf'})
  finally:engine.close()
 if folder and (folder/'正文/layout.json').exists():
  layout=read(folder/'正文/layout.json')
  if len(layout['pdf_info'])!=len(r.pages):raise ValueError('Layout page count differs from source PDF: '+str(pdf))
  for i,p in enumerate(layout['pdf_info']):
   for j,b in enumerate(p.get('preproc_blocks',[])):
    bi=b.get('index',j);kind=b.get('type','unknown');txt='\n'.join(x['content'] for x in walk(b) if isinstance(x.get('content'),str));imgs=list(dict.fromkeys(x['image_path'] for x in walk(b) if x.get('image_path')))
    asset=None
    if imgs:
     img=folder/'正文/images'/imgs[0]
     if img.exists():asset=add_asset(img,out,docid,f'block_p{i+1}_{bi}{img.suffix}')
    evidence.append({'evidence_id':f'{docid}:page{i+1}:block{bi}','document_id':docid,'page':i+1,'kind':kind,'text':txt or None,'asset':asset,'bbox':b.get('bbox'),'bbox_units':'pdf_points','provenance':'mineru_layout'})
  md=folder/'全文.md'
  if md.exists():
   title=re.search(r'^# (.+)$',md.read_text(),re.M);doc['title']=title.group(1) if title else None
 return doc

def empty_record(docid,rid,ev,kind='unresolved'):
 return {'record_id':f'{docid}:{rid}','document_id':docid,'source_reaction_id':rid,'record_kind':kind,'status':'candidate','completeness':'partial','evidence_ids':list(ev),'participants':[],'steps':[],'outcomes':[],'notes':[],'review':{'method':'none','reviewer':None,'checked_fields':[],'unresolved':['Experiment identity and chemical roles require source review.']}}

def import_chem(doc,folder,out):
 path=folder/'化学解析/demonstration_tables.json'
 if not path.exists():return []
 doc['chemical_raw_json']=add_asset(path,out,doc['document_id'],'mineru_raw_tables.json')
 raw=read(path);mols={m['mol_id']:m for m in raw['molecule_table']['data']};byev={e['evidence_id']:e for e in doc['evidence']};records=[]
 for rx in raw['reaction_table']['data']:
  page=rx.get('page_idx',0)+1;ev=f"{doc['document_id']}:page{page}:block{rx.get('block_index')}"
  if ev not in byev:ev=f"{doc['document_id']}:page{page}:text"
  rec=empty_record(doc['document_id'],rx['reaction_id'],[ev]);seen=set()
  for key,role in [('reactants_smiles','reactant'),('products_smiles','product')]:
   for obj in walk(rx.get(key,[])):
    if 'smiles' not in obj:continue
    mid=obj.get('mol_id');entity=(role,mid,obj.get('smiles'))
    if entity in seen:continue
    seen.add(entity);s=obj.get('smiles') or None
    rec['participants'].append({'entity_id':mid or f'unknown_{len(seen)}','role':role,'label':fact(mols.get(mid,{}).get('mol_identifier') or None,[ev]),'name':fact(None),'smiles':fact(s,[ev]),'smiles_parseable':smiles_ok(s),'structure_verified':False,'amounts':[]})
  for c in rx.get('reaction_conditions',[]):
   if c.get('type')=='molecule':
    mid=c.get('mol_id');m=mols.get(mid,{})
    rec['participants'].append({'entity_id':mid or f'condition_{len(rec["participants"])}','role':'unspecified_agent','label':fact(m.get('mol_identifier') or None,[ev]),'name':fact(None),'smiles':fact(m.get('mol_smiles') or None,[ev]),'smiles_parseable':smiles_ok(m.get('mol_smiles')),'structure_verified':False,'amounts':[]})
  # Retain unpaired OCR fragments; do not assign a reagent, temperature or stage by guesswork.
  rec['notes']=[fact(c.get('condition'),[ev]) for c in rx.get('reaction_conditions',[]) if c.get('type')=='text']
  rec['review']['unresolved']+=['Raw condition fragments are not assigned to experimental steps.','Missing yield/time in the API may be present in the source figure.','General scheme, prior-work example, mechanism and reported experiment are not yet distinguished.']
  if any(p['smiles_parseable'] is False for p in rec['participants']):rec['review']['unresolved'].append('Nonstandard/unexpanded or invalid SMILES; inspect aliases and source structure.')
  records.append(rec)
 return records

def prepare(args):
 out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=True)
 if (out/'dataset.json').exists():raise ValueError('Output already has a dataset; use another directory to preserve reviewed records.')
 dataset={'tool_version':VERSION,'schema_version':'1.0','created_at':datetime.now(timezone.utc).isoformat(),'documents':[],'candidates':[],'experiments':[],'run_log':[]};inputs=[]
 if args.manifest:
  for row in read(args.manifest):
   if args.chem_only and row.get('chem_status')!='已完成':continue
   inputs.append((Path(row['source']),Path(row['folder']),'main_text',None))
 elif args.folder:
  folder=Path(args.folder);pdfs=list(folder.glob('*.pdf'))
  if len(pdfs)!=1:raise ValueError('A document folder must contain exactly one top-level PDF.')
  inputs.append((pdfs[0],folder,'main_text',None))
 else:inputs.append((Path(args.pdf),None,'main_text',None))
 for pdf,folder,role,parent in inputs:
  doc=ingest(pdf,folder,out,role,parent,args.render)
  if any(d['document_id']==doc['document_id'] for d in dataset['documents']):
   dataset['run_log'].append({'event':'byte_identical_document_skipped','source_name':pdf.name,'document_id':doc['document_id']});continue
  dataset['documents'].append(doc)
  if folder:dataset['candidates'].extend(import_chem(doc,folder,out))
 if args.supplement:
  if len(dataset['documents'])!=1:raise ValueError('Explicit supplement association requires exactly one main document.')
  parent=dataset['documents'][0]['document_id'];dataset['documents'][0]['supplementary_material']={'status':'provided','note':'Explicitly associated by --supplement.'}
  for p in args.supplement:
   doc=ingest(Path(p),None,out,'supplementary',parent,args.render)
   if any(d['document_id']==doc['document_id'] for d in dataset['documents']):raise ValueError('Supplement duplicates an existing source PDF.')
   dataset['documents'].append(doc)
 write(out/'dataset.json',dataset)
 tasks=[]
 for d in dataset['documents']:
  for e in d['evidence']:
   # Include every figure/table and every page: fixed coverage, not question top-k only.
   if e['kind'] in ['page_text','image','table']:
    image_ev=next((v['evidence_id'] for v in d['evidence'] if v['page']==e['page'] and v['kind']=='page_image'),None)
    tasks.append({'task_id':e['evidence_id'],'document_id':d['document_id'],'evidence_id':e['evidence_id'],'visual_evidence_id':image_ev,'state':'pending','reason':'Identify experiments and cite only supplied evidence; page tasks also inspect prose and SI.'})
 write(out/'tasks.json',tasks);report(out)
 print(json.dumps({'documents':len(dataset['documents']),'candidates':len(dataset['candidates']),'tasks':len(tasks)},ensure_ascii=False))

def validate(data,base,strict_record_ids=()):
 from jsonschema import Draft202012Validator
 schema=read(ROOT/'references/experiment.schema.json');problems=[];docs={d['document_id']:d for d in data['documents']};ev={e['evidence_id']:e for d in docs.values() for e in d['evidence']};seen=set()
 if len(docs)!=len(data['documents']):problems.append('duplicate document_id')
 evidence_ids=[e['evidence_id'] for d in data['documents'] for e in d['evidence']]
 if len(evidence_ids)!=len(set(evidence_ids)):problems.append('duplicate evidence_id')
 for d in docs.values():
  pdf=base/d['source_pdf']
  if not pdf.resolve().is_relative_to(base.resolve()) or not pdf.exists() or sha(pdf)!=d['sha256']:problems.append(f"{d['document_id']}: source PDF hash mismatch")
  for e in d['evidence']:
   if e['document_id']!=d['document_id']:problems.append(e['evidence_id']+': evidence document mismatch')
   if not 1<=e['page']<=d['pages']:problems.append(e['evidence_id']+': page out of range')
   if e.get('asset') and (not (base/e['asset']).resolve().is_relative_to(base.resolve()) or not (base/e['asset']).is_file()):problems.append(e['evidence_id']+': asset missing')
 for rec in data['candidates']+data['experiments']:
  # A reviewed replacement may have the same id as its retained raw candidate.
  collection='experiments' if any(rec is x for x in data['experiments']) else 'candidates';key=(collection,rec.get('record_id'))
  if key in seen:problems.append('duplicate record id: '+str(key))
  seen.add(key)
  for err in Draft202012Validator(schema).iter_errors(rec):problems.append(f"{rec.get('record_id')}: {'/'.join(map(str,err.path))}: {err.message}")
  if not rec.get('record_id','').startswith(rec.get('document_id','')+':'):problems.append(str(key)+': record id must be namespaced by document id')
  if rec.get('document_id') not in docs:problems.append(str(key)+': unknown document')
  allowed={rec.get('document_id')}|{d['document_id'] for d in docs.values() if d.get('parent_document_id')==rec.get('document_id')}
  ids={p['entity_id'] for p in rec.get('participants',[])}
  if len(ids)!=len(rec.get('participants',[])):problems.append(str(key)+': duplicate participant entity_id')
  for s in rec.get('steps',[]):
   if any(x not in ids for x in s.get('participant_ids',[])):problems.append(str(key)+': step has unknown participant')
  for obj in walk(rec):
   for eid in obj.get('evidence_ids',[]):
    if eid not in ev or ev[eid]['document_id'] not in allowed:problems.append(str(key)+': invalid or unrelated evidence '+eid)
   if 'verification' in obj and 'value' in obj:
    if obj['value'] is not None and not obj.get('evidence_ids'):problems.append(str(key)+': non-null fact without evidence')
    if obj['value'] is None and obj['verification']!='missing':problems.append(str(key)+': null must be marked missing')
    if obj['value'] is not None and obj['verification']=='missing':problems.append(str(key)+': non-null fact cannot be missing')
    for error in quantity_errors(obj['value']):problems.append(str(key)+': '+error)
    if rec.get('status')=='reviewed' and obj.get('verification')=='unreviewed':problems.append(str(key)+': reviewed record contains unreviewed fact')
  if rec.get('status')=='reviewed' and rec.get('review',{}).get('method') not in ['human','assistant_visual']:problems.append(str(key)+': reviewed requires actual source review, not model consensus')
  for p in rec.get('participants',[]):
   if p.get('structure_verified') and p['smiles']['value'] is not None and smiles_ok(p['smiles']['value']) is False:problems.append(str(key)+': verified structure has invalid standard SMILES')
  for error in record_errors(rec,strict=rec['record_id'] in strict_record_ids):problems.append(str(key)+': '+error)
  if rec.get('completeness')=='complete':
   supplied_si=any(d.get('role')=='supplementary' for d in docs.values() if d.get('parent_document_id')==rec.get('document_id'))
   si=docs.get(rec.get('document_id'),{}).get('supplementary_material',{})
   if not supplied_si and not (si.get('status')=='not_applicable' and si.get('note')):problems.append(str(key)+': complete record requires SI or explicit not_applicable justification')
 return sorted(set(problems))

def merge(args):
 base=Path(args.dataset).resolve().parent;data=read(args.dataset);new=read(args.records)
 if not isinstance(new,list):raise ValueError('Records file must be a JSON array.')
 if len({r['record_id'] for r in new})!=len(new):raise ValueError('Duplicate submitted record_id')
 trial=copy.deepcopy(data);indexed={r['record_id']:r for r in trial['experiments']}
 for r in copy.deepcopy(new):
  if indexed.get(r['record_id'],{}).get('status')=='reviewed':
   raise ValueError('Cannot overwrite a source-reviewed record through machine merge')
  r['status']='needs_review'
  r['review']['method']='model'
  r['review']['checked_fields']=[]
  for obj in walk(r):
   if 'verification' in obj and 'value' in obj:
    obj['verification']='missing' if obj['value'] is None else 'unreviewed'
  for participant in r.get('participants',[]):
   participant['structure_verified']=False
   value=participant['smiles']['value']
   participant['smiles_parseable']=smiles_ok(value) if isinstance(value,str) else None
  normalize_record(r)
  indexed[r['record_id']]=r
 trial['experiments']=list(indexed.values());errors=validate(trial,base,strict_record_ids={r['record_id'] for r in new})
 if errors:raise ValueError('Merge refused: '+json.dumps(errors,ensure_ascii=False))
 trial['tool_version']=VERSION
 trial['run_log'].append({'event':'records_merged','count':len(new),'at':datetime.now(timezone.utc).isoformat()});write(args.dataset,trial);report(base);print('Merged',len(new),'records')

def report(base):
 d=read(base/'dataset.json');ev={x['evidence_id']:x for doc in d['documents'] for x in doc['evidence']};counts={'documents':len(d['documents']),'candidates':len(d['candidates']),'structured_records':len(d['experiments']),'source_reviewed_records':sum(r['status']=='reviewed' for r in d['experiments']),'semantic_coverage_complete':False}
 write(base/'summary.json',counts);lines=['# 化学反应提取结果','',f"文献{counts['documents']}份；原始候选{counts['candidates']}条；结构化记录{counts['structured_records']}条；原文核查记录{counts['source_reviewed_records']}条。",'','页码和图表已枚举不代表全部反应已提取。未核查候选不得作为可靠实验事实。','']
 for r in d['experiments']:
  lines += ['## '+r['record_id'],'',f"类型：{r['record_kind']}；状态：{r['status']}；完整性：{r['completeness']}",'']
  labels={p['entity_id']:str(p['label']['value'] or p['name']['value'] or p['entity_id']) for p in r['participants']}
  def display_field(field):
   fact=field['fact'];targets=field.get('target_entity_ids',[])
   target='；关联对象：'+', '.join(labels.get(t,t)+' ['+t+']' for t in targets) if targets else ''
   refs='；字段依据：'+', '.join(fact.get('evidence_ids',[])) if fact.get('evidence_ids') else ''
   return field['name']+': '+json.dumps(fact['value'],ensure_ascii=False)+target+'；核查：'+fact.get('verification','missing')+refs
  for p in r['participants']:
   lines.append(f"- {p['role']}: {labels[p['entity_id']]} [{p['entity_id']}]; SMILES候选={p['smiles']['value'] or '未确认'}；结构核实={p.get('structure_verified',False)}")
   for a in p['amounts']:lines.append('  - '+display_field(a))
  for s in r['steps']:
   lines.append('\n步骤 '+s['step_id'])
   if s.get('participant_ids'):lines.append('涉及对象：'+', '.join(labels.get(t,t)+' ['+t+']' for t in s['participant_ids']))
   for c in s['conditions']:lines.append('- '+display_field(c))
  for o in r['outcomes']:lines.append('- '+display_field(o))
  for note in r['notes']:lines.append('\n备注：'+str(note['value']))
  for eid in r['evidence_ids']:
   e=ev[eid];lines.append(f"\n依据：PDF第{e['page']}页，{eid}"+(f"；[原图]({e['asset']})" if e.get('asset') else ''))
  lines.append('\n待核查：'+'；'.join(r['review']['unresolved']))
 (base/'文献实验记录.md').write_text('\n'.join(lines)+'\n')
 db=sqlite3.connect(base/'experiments.sqlite');db.execute('CREATE TABLE IF NOT EXISTS experiments (record_id TEXT PRIMARY KEY, document_id TEXT, status TEXT, content TEXT)');db.execute('DELETE FROM experiments')
 for r in d['experiments']:db.execute('INSERT INTO experiments VALUES (?,?,?,?)',(r['record_id'],r['document_id'],r['status'],json.dumps(r,ensure_ascii=False)))
 db.commit();db.close()

def main():
 if len(sys.argv)>1 and sys.argv[1]=='standard':
  from standard_records import main as standard_main
  return standard_main(sys.argv[2:])
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
 a=sub.add_parser('prepare');g=a.add_mutually_exclusive_group(required=True);g.add_argument('--manifest');g.add_argument('--folder');g.add_argument('--pdf');a.add_argument('--supplement',action='append');a.add_argument('--chem-only',action='store_true');a.add_argument('--render',action='store_true');a.add_argument('--out',required=True);a.set_defaults(func=prepare)
 a=sub.add_parser('merge');a.add_argument('--dataset',required=True);a.add_argument('--records',required=True);a.set_defaults(func=merge)
 a=sub.add_parser('validate');a.add_argument('--dataset',required=True);a.add_argument('--strict',action='store_true');a.set_defaults(func=lambda a:check(a))
 a=sub.add_parser('audit');a.add_argument('--dataset',required=True);a.set_defaults(func=audit)
 a=sub.add_parser('crop');a.add_argument('--dataset',required=True);a.add_argument('--evidence',required=True);a.add_argument('--box',nargs=4,type=int,required=True);a.add_argument('--out',required=True);a.set_defaults(func=crop)
 a=sub.add_parser('inspect');a.add_argument('--dataset',required=True);a.add_argument('--document',required=True);a.add_argument('--page',type=int,required=True);a.set_defaults(func=inspect)
 a=sub.add_parser('query');a.add_argument('--dataset',required=True);a.add_argument('--text',default='');a.add_argument('--include-unreviewed',action='store_true');a.add_argument('--category',choices=['yield','recovery','conversion','time','temperature','percentage','mass','volume','equivalents','amount','concentration','selectivity','measurement_basis','atmosphere','workup','purification','operation','solvent','other']);a.set_defaults(func=query)
 args=p.parse_args()
 if hasattr(args,'limit') and args.limit<0:p.error('--limit must be nonnegative (0 = all remaining tasks)')
 try:args.func(args)
 except Exception as e:print('ERROR:',str(e),file=sys.stderr);return 1
 return 0

def check(a):
 data=read(a.dataset);base=Path(a.dataset).resolve().parent
 errors=validate(data,base,strict_record_ids={r['record_id'] for r in data['experiments']} if a.strict else ())
 if 'standard' in data:
  from standard_records import validate as standard_validate
  errors += standard_validate(data['standard'],base)
 print(json.dumps({'valid':not errors,'errors':errors,'standard_v3_present':'standard' in data,'chemical_accuracy_verified':False},ensure_ascii=False,indent=2))
 if errors:raise ValueError('Validation failed')
def audit(a):
 data=read(a.dataset)
 result=dataset_audit(data)
 print(json.dumps(result,ensure_ascii=False,indent=2))

def dataset_audit(data):
 from collections import Counter
 return {'tool_version':VERSION,'chemical_accuracy_verified':False,
         'record_kind_counts':dict(Counter(r['record_kind'] for r in data['experiments'])),
         'source_reviewed_records':sum(r.get('status')=='reviewed' for r in data['experiments']),
         'documents_without_linked_si':[d['document_id'] for d in data['documents'] if d.get('role')!='supplementary' and not any(s.get('role')=='supplementary' and s.get('parent_document_id')==d['document_id'] for s in data['documents'])],
         'limitations':'Advisory structural checks only. Missing SI, structures or methods may reflect the source. No warning does not establish correctness, recall or suitability for expert rules.',
         'warnings':[{'record_id':r['record_id'],'record_kind':r['record_kind'],'evidence_ids':r.get('evidence_ids',[]),'warnings':w} for r in data['experiments'] if (w:=warnings(r))]}

def inspect(a):
 d=read(a.dataset);doc=next((x for x in d['documents'] if x['document_id']==a.document),None)
 if doc is None or not 1<=a.page<=doc['pages']:raise ValueError('Unknown document or page out of range')
 evidence=[e for e in doc['evidence'] if e['page']==a.page];ids={e['evidence_id'] for e in evidence}
 print(json.dumps({'document_id':doc['document_id'],'source_pdf':doc['source_pdf'],'page':a.page,'evidence':evidence,'candidate_hints':[r for r in d['candidates'] if ids.intersection(r['evidence_ids'])]},ensure_ascii=False,indent=2))

def query(a):
 from record_quality import fields
 d=read(a.dataset);selected=[r for r in d['experiments'] if (a.include_unreviewed or r['status']=='reviewed') and a.text.casefold() in json.dumps(r,ensure_ascii=False).casefold()]
 if getattr(a,'category',None):
  selected=[{'record':r,'matched_fields':matched} for r in selected if (matched:=[dict(f,category=category(f['name']),stored_category=f.get('category')) for f in fields(r) if category(f['name'])==a.category])]
 print(json.dumps(selected,ensure_ascii=False,indent=2))
def crop(a):
 from PIL import Image
 base=Path(a.dataset).resolve().parent
 data=read(a.dataset)
 e=next((e for d in data['documents'] for e in d['evidence'] if e['evidence_id']==a.evidence),None)
 if not e or not e.get('asset'):raise ValueError('Evidence must identify an image asset')
 asset=(base/e['asset']).resolve()
 if not asset.is_relative_to(base):raise ValueError('Image outside dataset')
 with Image.open(asset) as im:
  x0,y0,x1,y1=a.box
  if not 0<=x0<x1<=im.width or not 0<=y0<y1<=im.height:raise ValueError('Crop outside image; box uses image pixels')
  dest=Path(a.out);dest.parent.mkdir(parents=True,exist_ok=True);im.crop(a.box).save(dest)
 write(str(dest)+'.evidence.json',{'evidence_id':a.evidence,'source_asset':e['asset'],'source_sha256':sha(asset),'page':e['page'],'box_pixels':a.box,'crop_sha256':sha(dest)})
 print(dest)

if __name__=='__main__':sys.exit(main())
