"""Local, row-level comparison of an independent report with saved Promethease.

Never imported by inference. All outputs must be outside the release folder.
Aliases are established from the public catalogue, not from co-located probes.
"""
import argparse
import base64
from collections import Counter, defaultdict
import hashlib
import html
import json
from pathlib import Path
import re
import sqlite3
import zlib


def norm(value):
    return ' '.join(html.unescape(re.sub('<[^>]*>', '',str(value or ''))).split()).casefold()


def genotype(value):
    return tuple(sorted(str(value).upper().strip('()').split(';')))


def read_rows(path):
    chunks=re.findall(r"mygenos\.push\.apply\(mygenos,decompressString\('([^']+)'\)\)",Path(path).read_text())
    if not chunks:raise ValueError('No saved Promethease records found.')
    return [r for c in chunks for r in json.loads(zlib.decompress(base64.b64decode(c)))]


def compare_rows(rows,report,public_aliases):
    k=report['knowledge'];findings={r['id']:r for r in k['findings']}
    unresolved={r['id']:r['reason'] for r in k['unresolved']}
    rule_audit={r['id']:r for r in k['genoset_audit']}
    out=[];counts=Counter();reasons=Counter();content=Counter();strata=defaultdict(Counter)
    for index,row in enumerate(rows):
        key=row.get('rsnum',row.get('title','')).lower();variant='rsnum' in row
        result=findings.get(key);alias=None;reason=''
        if variant:
            if result:status='exact_genotype' if genotype(result.get('source_genotype'))==genotype(row.get('geno')) else 'genotype_disagreement'
            elif public_aliases.get(key) in findings:
                alias=public_aliases[key];result=findings[alias];status='explicit_alias_represented'
            else:status='unresolved';reason=unresolved.get(key,unresolved.get(public_aliases.get(key),'not_resolved_from_raw_input'))
            if status=='genotype_disagreement' and 'literal_source_indel_code_match' in result.get('flags',[]):
                reason='indel_code_vs_sequence_notation_not_reconciled'
            if status=='genotype_disagreement' and result.get('sequence_resolution'):
                resolved=result['sequence_resolution'];baseline=genotype(row.get('geno'))
                if len(baseline)==2 and set(baseline)<=set('ID'):
                    mapped=tuple(sorted(resolved['definition']['source_insertion'] if a=='I' else '-' for a in baseline))
                    if mapped==genotype(result['source_genotype']):
                        status='public_sequence_notation_equivalent';reason='Explicit public sequence mapping; not literal genotype equality or clinical validation.'
            if status=='genotype_disagreement' and result.get('curated_sequence_evidence'):
                definition=result['curated_sequence_evidence'];raw=result.get('raw_genotype','');orientation=definition.get('orientation')
                translated=raw.translate(str.maketrans('ACGT','TGCA')) if orientation=='minus' else raw
                if len(raw)==1:raw*=2;translated*=2
                if orientation in {'plus','minus'} and tuple(sorted(raw))==genotype(row.get('geno')) and tuple(sorted(translated))==genotype(result.get('source_genotype')):
                    status='public_strand_notation_equivalent';reason='The historical row uses the raw forward genotype; the source-page genotype uses the independently corroborated strand. Not literal genotype equality or clinical validation.'
        else:
            status='matched_rule' if result else 'unresolved_rule'
            reason=rule_audit.get(key,{}).get('reason','rule_absent') if not result else ''
        counts[status]+=1
        if reason:reasons[reason]+=1
        mag=float(row.get('magnitude') or 0);band='magnitude >= 3' if mag>=3 else 'magnitude 1 to <3' if mag>=1 else 'magnitude <1'
        strata[band][status]+=1
        item=dict(row=index,id=key,kind='variant' if variant else 'genoset',status=status,reason=reason,alias_target=alias,
                  baseline_genotype=row.get('geno'),app_genotype=result.get('source_genotype') if result else None,
                  baseline_summary=row.get('genosummary',''),app_summary=result.get('summary','') if result else '',
                  baseline_magnitude=mag,app_magnitude=result.get('magnitude') if result else None,
                  baseline_repute=row.get('repute','Unknown'),app_repute=result.get('repute') if result else None)
        if status in {'exact_genotype','matched_rule'}:
            for field,a,b in [('summary',row.get('genosummary'),result.get('summary') if result.get('has_genotype_summary') else ''),('repute',row.get('repute','Unknown'),result.get('repute','Unknown')),('magnitude',str(mag),str(float(result.get('magnitude') or 0)))]:
                equal=norm(a)==norm(b);item[field+'_equal']=equal;content[field+('_same' if equal else '_different')]+=1
            if row.get('genosummary'):
                content['baseline_explained_rows_matched']+=1
                content['baseline_explanations_same' if item['summary_equal'] else 'baseline_explanations_changed']+=1
                if not item['summary_equal'] and norm(row['genosummary']).replace('"',"'")==norm(result['summary']).replace('"',"'"):
                    content['changed_explanations_quotes_only']+=1
            if mag==1 and result.get('magnitude',0)==0:content['magnitude_1_vs_0_default_difference']+=1
            if norm(row.get('repute','Unknown'))=='unknown' and not result.get('repute'):content['unknown_vs_empty_label']+=1
        if status=='explicit_alias_represented' and row.get('genosummary'):
            content['alias_baseline_explained_rows']+=1
            content['alias_summary_same' if norm(row['genosummary'])==norm(result['summary']) else 'alias_summary_different']+=1
        if not variant:item['rule_audit']=rule_audit.get(key,{})
        out.append(item)
    represented={r['id'] for r in out if r['status'] in {'exact_genotype','matched_rule','public_sequence_notation_equivalent','public_strand_notation_equivalent'}}|{r['alias_target'] for r in out if r['alias_target']}
    extras=[dict(id=r['id'],kind=r['kind'],summary=r['summary'],magnitude=r['magnitude']) for r in k['findings'] if r['id'] not in represented]
    return dict(baseline_rows=len(rows),baseline_variant_rows=sum('rsnum' in r for r in rows),baseline_explained_rows=sum(bool(r.get('genosummary')) for r in rows),
                counts=dict(counts),content_agreement=dict(content),unresolved_reasons=dict(reasons),magnitude_strata={k:dict(v) for k,v in strata.items()},
                app_counts=k['counts'],app_additional_records=len(extras),extra_genosets=[r for r in extras if r['kind']=='genoset'],rows=out,
                limits=['One person and one historical report; not a current-service or clinical-accuracy validation.',
                        'Exact matches require the same marker ID and genotype. Public aliases are counted separately; their interpretation equivalence is not assumed.',
                        'Text agreement ignores HTML tags, case and whitespace only. Changed text is not automatically an error.',
                        'Extra records are not proof of additional useful discoveries. Source magnitude is not clinical risk.'])


def main():
    p=argparse.ArgumentParser();p.add_argument('--promethease',type=Path,required=True);p.add_argument('--report',type=Path,required=True);p.add_argument('--knowledge',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.resolve().is_relative_to(Path(__file__).resolve().parents[1]):p.error('Keep personal comparisons outside the release folder.')
    db=sqlite3.connect(a.knowledge.resolve().as_uri()+'?mode=ro',uri=True);aliases={key:r['alias'] for key,raw in db.execute('select id,data from markers') if (r:=json.loads(raw)).get('alias')};db.close()
    d=compare_rows(read_rows(a.promethease),json.loads(a.report.read_text()),aliases)
    d['input_sha256']={label:hashlib.sha256(path.read_bytes()).hexdigest() for label,path in [('baseline',a.promethease),('standalone_report',a.report),('public_knowledge',a.knowledge)]}
    a.out.mkdir(parents=True,exist_ok=True);(a.out/'comparison.json').write_text(json.dumps(d,indent=2))
    print(json.dumps({k:v for k,v in d.items() if k not in {'rows'}},indent=2))

if __name__=='__main__':main()
