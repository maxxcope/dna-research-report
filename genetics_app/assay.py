"""Preserve attributed assay cautions without inheriting genotype conclusions."""
import re

CAUTION=re.compile(r'miscall|genotyping error|unreliable (?:probe|genotyp|call)|probe[^.\n]{0,100}\bunreliable\b|error[- ]prone|prone to false results|false[- ]positive (?:genotyp|call)',re.I)


def warning_passage(text):
    match=CAUTION.search(text or '')
    if not match:return ''
    start=max(text.rfind('\n',0,match.start()),text.rfind('. ',0,match.start())+1,0)
    end=text.find('\n',match.end())
    if end<0:end=len(text)
    return text[start:min(end,match.end()+500)].strip()


def probe_warning_records(markers,genotypes):
    """Index warnings on public probe pages, regardless of genotype matching.

    This is cautionary context only. No source genotype, disease assertion,
    magnitude or clinical classification is inherited from a probe page.
    """
    out={}
    for key,record in list(markers.items())+list(genotypes.items()):
        probe=key.split('(')[0]
        if not probe.startswith('i'):continue
        text=record.get('text','')
        quote=warning_passage(text)
        if not quote and record.get('summary'):quote=warning_passage(record['summary'])
        if not quote:continue
        # A page may continue with a conditional disease claim after a warning.
        # Retain only the warning clause as probe-level context.
        match=CAUTION.search(quote)
        boundaries=list(re.finditer(r';|,?\s+but\b',quote[:match.start()],re.I))
        if boundaries:quote=quote[boundaries[-1].end():].strip()
        match=CAUTION.search(quote);end=re.search(r';|,?\s+but\b|\s+would otherwise\b',quote[match.end():],re.I)
        if end:quote=quote[:match.end()+end.start()].strip()
        item=dict(probe_id=probe,source_id=key,source=record['source'],warning=quote,
            scope='Public warning about this input probe or some of its source genotype pages; the page genotype is not asserted to match this call.')
        if item not in out.setdefault(probe,[]):out[probe].append(item)
    return out


def finding_warning_records(finding):
    rows=[]
    for key in ('summary','genotype_text','variant_text'):
        quote=warning_passage(finding.get(key,''))
        if quote:
            named_probes=set(re.findall(r'\bi\d+\b',quote,re.I))
            if named_probes and 'raw_ids' in finding and not named_probes&set(finding['raw_ids']):continue
            source=finding.get('source','')
            if key=='variant_text':source=source.split('(',1)[0]
            rows.append(dict(source_id=finding['id'],source=source,warning=quote,
                scope='Warning in the matched source record or its general variant context.'))
            break
    rows.extend(finding.get('probe_assay_warnings',[]))
    return list({(r['source_id'],r['warning']):r for r in rows}.values())
