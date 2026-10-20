"""Finite evidence selections from a saved report, without LLMs or network calls."""
from collections import defaultdict, deque
import re

from .navigation import COLLECTIONS, MARKER_TEXT

SUPPORT = {'/knowledge/findings', '/knowledge/unresolved', '/knowledge/genoset_audit',
           '/clinvar', '/traits', '/exclusions', '/finding_audit/flagged_findings'}
WARNING = {'/knowledge/unresolved', '/knowledge/genoset_audit', '/exclusions',
           '/finding_audit/flagged_findings'}


def strings(value):
    if isinstance(value, str):
        return [value]
    if isinstance(value, list):
        return [s for s in value if isinstance(s, str)]
    return []


def primary_ids(value):
    result = set()
    for key in ('id', 'rsid', 'rs', 'marker_ids'):
        for text in strings(value.get(key)):
            if re.fullmatch(r'(?:rs|gs|i)\d+', text):
                result.add(text)
    return result


def required_ids(value):
    """Forward references only; shared IDs never make another rule required."""
    result = set()
    identity = value.get('id', '')
    if isinstance(identity, str) and identity.startswith(('conditional:', 'assay:')):
        candidate = identity.split(':', 1)[1]
        if re.fullmatch(r'(?:rs|gs|i)\d+', candidate):
            result.add(candidate)
    for field in ('dependencies', 'missing_markers', 'blocking_markers'):
        for text in strings(value.get(field)):
            result.update(MARKER_TEXT.findall(text))
    if isinstance(value.get('rule'), str):
        result.update(MARKER_TEXT.findall(value['rule']))
    for field in ('variants', 'source_patterns', 'observed_support', 'unobserved_tests'):
        children = value.get(field, [])
        if isinstance(children, dict):
            result.update(key for key in children if re.fullmatch(r'(?:rs|gs|i)\d+', key))
            children = [children]
        for child in children if isinstance(children, list) else []:
            if isinstance(child, dict):
                result.update(primary_ids(child))
                result.update(strings(child.get('marker')))
    observed = value.get('observed', {})
    if isinstance(observed, dict):
        result.update(strings(observed.get('marker')))
    return {x for x in result if re.fullmatch(r'(?:rs|gs|i)\d+', x)}


def probe_ids(value):
    result = set()
    for field in ('raw_ids', 'public_aliases', 'colocated_raw_ids', 'marker_ids', 'input_markers'):
        result.update(strings(value.get(field)))
    calls = value.get('dependency_calls', {})
    if isinstance(calls, dict):
        for call in calls.values():
            if isinstance(call, dict):
                result.update(probe_ids(call))
    return result


def gene_ids(value):
    """Exact declared gene identifiers, not names guessed from source prose."""
    result = set()
    for field in ('gene', 'reported_gene'):
        for text in strings(value.get(field)):
            for gene in re.split(r'[|,;/\s]+', text):
                gene = re.sub(r':\d+$', '', gene)
                if re.fullmatch(r'[A-Za-z][A-Za-z0-9_.-]*', gene):
                    result.add(gene.upper())
    return result


class EvidenceIndex:
    """Search metadata, then select exact source records with forward support."""
    def __init__(self, report):
        from .bundle import at_path
        self.report = report
        self.records = {}
        self.by_id = defaultdict(set)
        self.by_gene = defaultdict(set)
        self.warnings = defaultdict(set)
        self.collection = {}
        for prefix in COLLECTIONS:
            try:
                items = at_path(report, prefix)
            except (KeyError, IndexError, TypeError):
                continue
            if not isinstance(items, list):
                continue
            for i, value in enumerate(items):
                path = f'{prefix}/{i}'
                self.records[path] = value
                self.collection[path] = prefix
                if not isinstance(value, dict):
                    continue
                if prefix in SUPPORT:
                    for marker in primary_ids(value):
                        self.by_id[marker].add(path)
                    for gene in gene_ids(value):
                        self.by_gene[gene].add(path)
                if prefix in WARNING:
                    for marker in primary_ids(value) | probe_ids(value):
                        self.warnings[marker].add(path)

    def search(self, query):
        """Return exact-word metadata matches, with existing summaries first."""
        query = query.strip()
        if not query:
            raise ValueError('Enter a gene, marker or source-label phrase.')
        pattern = re.compile(r'(?<![\w])' + re.escape(query) + r'(?![\w])', re.I)
        results = []
        for path, value in self.records.items():
            if not isinstance(value, dict):
                continue
            text = ' '.join(str(value.get(k, '')) for k in ('id', 'rsid', 'rs', 'gene', 'reported_gene', 'title', 'summary', 'conditions', 'trait'))
            if pattern.search(text):
                label = value.get('title') or value.get('summary') or value.get('conditions') or value.get('id') or value.get('rsid') or path
                results.append({'path': path, 'label': label})
        return sorted(results, key=lambda r: (not r['path'].startswith('/conclusions/cards/'), r['path']))

    def select(self, paths):
        """Closed required set. No reverse rule links and no co-packed neighbors."""
        seeds = list(dict.fromkeys(paths))
        if not seeds or any(path not in self.records for path in seeds):
            raise ValueError('Select existing record paths from the report search results.')
        pending = deque(seeds)
        selected = set()
        unresolved = set()
        while pending:
            path = pending.popleft()
            if path in selected:
                continue
            selected.add(path)
            value = self.records[path]
            if not isinstance(value, dict):
                continue
            # A gene-labelled medication summary must not hide other measured
            # alleles in that gene just because they are outside its matched rule.
            identity = value.get('id', '')
            if path.startswith('/conclusions/cards/') and isinstance(identity, str) and identity.startswith('medication:'):
                gene = identity.split(':', 1)[1].upper()
                pending.extend(sorted(self.by_gene.get(gene, set()) - selected))
            for marker in primary_ids(value) | required_ids(value):
                matches = self.by_id.get(marker, set())
                if not matches:
                    unresolved.add(marker)
                pending.extend(sorted(matches - selected))
            for probe in primary_ids(value) | probe_ids(value):
                pending.extend(sorted(self.warnings.get(probe, set()) - selected))
        from .bundle import at_path
        values = {path: self.records[path] for path in sorted(selected)}
        for path in ('/limits', '/conclusions/evidence_limits', '/conclusions/method'):
            try:
                values[path] = at_path(self.report, path)
            except (KeyError, IndexError, TypeError):
                pass
        return {'scope': seeds, 'records': values,
                'no_separate_record_for': sorted(unresolved),
                'selection_rule': 'Selected records, their forward marker/rule references, saved clinical context and linked warning records. Gene-labelled medication summaries also include directly annotated records for that gene, including unresolved calls. No reverse expansion to other rules or shared-file neighbors.',
                'stop': 'This is the complete saved-record set for this selection rule. Other subjects are optional. This is not all biological knowledge about a gene; absence of a separate record is not absence of a genotype or condition.'}


def write_focus_views(report, folder, maps, part_bytes, qualifications):
    """Prebuild finite reading sets for existing evidence summaries, not new claims."""
    import hashlib
    import os
    from pathlib import Path
    from .bundle import compact, compact_part, expand_part, records, render_part, write_text
    index = EvidenceIndex(report)
    views = {}
    inventory = []
    specs = []
    for seed, value in index.records.items():
        if seed.startswith('/conclusions/cards/'):
            selected = index.select([seed])
            specs.append(dict(selected, seed=seed, label=value.get('title', seed), area=None, details=[]))
    from .research_navigation import selections
    specs.extend(selections(report))
    for selection in specs:
        seed = selection['seed']
        is_medication = seed in index.records and str(index.records[seed].get('id', '')).startswith('medication:')
        digest = hashlib.sha256(seed.encode()).hexdigest()[:20]
        map_path = f'maps/focus/{digest}.md'
        rows = []
        for path, source in selection['records'].items():
            limit = min(8000, (part_bytes - 8000) // 2)
            if len(compact(source)) <= limit:
                rows.append({'path': path, 'kind': 'value', 'value': source})
                continue
            for _, row in records(source, record_bytes=limit):
                row = dict(row, path=path + row['path'])
                if row['kind'] == 'alias':
                    row['target'] = path + row['target']
                rows.append(row)
        files = []

        def save(batch):
            if not batch:
                return
            payload = compact_part(batch)
            if compact(expand_part(payload)) != compact(batch):
                raise ValueError('Focused export failed its lossless check.')
            name = f'focus/{digest}-{len(files) + 1:04d}.md'
            relative = os.path.relpath(map_path, Path(name).parent).replace(os.sep, '/')
            from .navigation import safe_label
            text = render_part('findings', payload, qualifications)
            text = text.replace('# Related findings and warnings', '# ' + safe_label(selection['label'], 180), 1)
            text = text.replace('Use index.tsv to locate related paths before interpreting fragments.',
                                'All fragments and alias targets for this selection are in its required reading files; collect those files before interpreting them.')
            text = text.replace('## Reading the data',
                                f'## Selected evidence only\n\n'
                                f'Read every file in the [required reading plan]({relative}). '
                                'Then stop: this selection has no additional required companion branches. '
                                'Other topics in the full report are optional. Preserve the selected result’s status and limitations.\n\n'
                                '## Reading the data', 1)
            if len(text.encode()) > part_bytes:
                if len(batch) == 1:
                    raise ValueError('A focused record exceeds the reading-file budget.')
                middle = len(batch) // 2
                save(batch[:middle]); save(batch[middle:])
                return
            write_text(Path(folder) / name, text)
            files.append(name)
            inventory.append({'file': name, 'scope': seed, 'bytes': len(text.encode()),
                              'sha256': hashlib.sha256(text.encode()).hexdigest()})
        batch = []
        weight = 0
        for row in rows:
            size = len(compact(row))
            if batch and weight + size > part_bytes - 8000:
                save(batch);batch=[];weight=0
            batch.append(row);weight+=size
        save(batch)
        label = selection['label']
        detail = ('Medication summaries include other directly annotated records for that gene. '
                  'Different source labels may describe individual markers or partial patterns; do not combine them into a complete medication phenotype. ') if is_medication else ''
        maps.tree(map_path, label, [(f'Required section {i + 1} of {len(files)}', name) for i, name in enumerate(files)],
                  'Read all listed sections for this selected result, including its supporting records, status, warnings and limitations. '
                  'Then stop; there are no required branches beyond this plan. Other subjects and detailed-data routes are optional. '
                  + detail +
                  'Map labels do not add interpretation or change the original result. '
                  'The complete original report remains in full-report.json.gz.')
        views[seed] = {'map': map_path, 'files': files, 'record_paths': list(selection['records']),
                       'selection_rule': selection['selection_rule'],
                       'no_separate_record_for': selection.get('no_separate_record_for', []),
                       'label': label, 'area': selection['area'], 'details': selection['details']}
    return views, inventory
