#!/usr/bin/env python3
"""Read-only corpus triage; write a separate review queue, never repair facts."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

from reaction_tools import dataset_audit


def scan(root, excluded=()):
    root = Path(root).resolve()
    seen = {}
    result = {'chemical_accuracy_verified': False, 'datasets': [], 'duplicates': [], 'errors': []}
    kinds, issues, versions = Counter(), Counter(), Counter()
    for path in sorted(root.rglob('dataset.json')):
        if any(part in excluded for part in path.relative_to(root).parts):
            continue
        try:
            raw = path.read_bytes()
            data = json.loads(raw)
            # Same document can have several folder mappings. Never silently deduplicate conflicting outputs.
            identity = tuple(sorted(d['document_id'] for d in data['documents']))
            digest = hashlib.sha256(raw).hexdigest()
            key = (identity, digest)
            if key in seen:
                result['duplicates'].append({'path': str(path), 'same_as': seen[key]})
                continue
            if any(k[0] == identity for k in seen):
                result['errors'].append({'path': str(path), 'error': 'same document has differing dataset versions; both retained for review'})
            seen[key] = str(path)
            audit = dataset_audit(data)
            audit.update(dataset=str(path), dataset_sha256=digest,
                         extraction_version=data.get('extraction_profile', {}).get('version', 'legacy'))
            result['datasets'].append(audit)
            kinds.update(audit['record_kind_counts'])
            versions[audit['extraction_version']] += 1
            for record in audit['warnings']:
                issues.update(set(w.split(':', 1)[0] for w in record['warnings']))
        except (ValueError, KeyError, TypeError, OSError) as error:
            result['errors'].append({'path': str(path), 'error': str(error)})
    result['summary'] = {'dataset_versions': dict(versions), 'record_kind_counts': dict(kinds),
                         'flagged_records_by_issue': dict(issues),
                         'datasets_scanned': len(result['datasets']),
                         'records_scanned': sum(kinds.values()),
                         'source_reviewed_records': sum(d['source_reviewed_records'] for d in result['datasets']),
                         'excluded_directory_names': list(excluded)}
    return result


def save(result, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    (out/'audit.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    with (out/'复核队列.tsv').open('w', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['数据集', '提取版本', '记录ID', '类型', '问题', '证据ID'])
        for dataset in result['datasets']:
            for row in dataset['warnings']:
                writer.writerow([dataset['dataset'], dataset['extraction_version'], row['record_id'],
                                 row['record_kind'], ' | '.join(row['warnings']), ', '.join(row['evidence_ids'])])
    lines = ['# 机器记录复核线索', '', '本清单不修改任何原始记录。警告是复核线索，不是错误率；无警告不代表化学正确或提取完整。',
             '同文献同内容的重复目录去重；同文献不同内容的版本均保留并列入冲突。人工修订目录按调用参数排除。', '',
             '扫描数据集：'+str(result['summary']['datasets_scanned']),
             '记录总数（含非实验）：'+str(result['summary']['records_scanned']), '', '## 按记录类型', '']
    lines += [f'- {k}: {v}' for k, v in result['summary']['record_kind_counts'].items()]
    lines += ['', '## 待核查线索（类别可重叠）', '']
    lines += [f'- {k}: {v}' for k, v in result['summary']['flagged_records_by_issue'].items()]
    lines += ['', '逐条位置、记录ID与证据见复核队列.tsv；完整扫描信息与原数据哈希见audit.json。',
              '优先核对 characterization_only、reaction_identity、outcome_target；再核对测定口径与结构缺项。',
              '对整张原图/表另做遗漏检查，本工具不能从已有JSON证明反应召回率。',
              '无反应实验不强制产物；缺失条件保留未知；未关联SI不等于SI不存在。']
    (out/'复核摘要.md').write_text('\n'.join(lines)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--exclude', action='append', default=[])
    args = parser.parse_args()
    root, out = Path(args.root).resolve(), Path(args.out).resolve()
    if not root.is_dir():
        parser.error('root must be an existing directory')
    if out == root or root.is_relative_to(out):
        parser.error('output must not contain the source corpus')
    if out.exists():
        parser.error('output must be a new directory')
    result = scan(root, args.exclude)
    save(result, out)
    print(json.dumps(result['summary'], ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
