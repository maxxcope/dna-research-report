"""Readable research routes, retaining bulk reference data as optional detail."""
import re


TITLES = {
    'array_quality': ('methods', 'Array quality and flagged positions'),
    'frequencies': ('population', 'Reference genotype frequencies'),
    'paternal_lineage': ('lineages', 'Paternal lineage candidates and limitations'),
    'maternal_lineage': ('lineages', 'Maternal lineage candidates and limitations'),
    'score_coverage': ('scores', 'Measured score coverage'),
    'imputed_score': ('scores', 'Estimated score, components and validation'),
    'phase': ('phase', 'Targeted phase, observations and limitations'),
}


def selections(report):
    """Select original paths, never recompute or abbreviate result values."""
    from .bundle import child_path
    research = report.get('research', {})
    if not isinstance(research, dict):
        return []
    result = []

    def add(seed, value, area, label, skip=()):
        # Unknown fields are included by default. Only known bulk fields move
        # to detail routes; failed/unsupported status and explanations stay.
        values = {child_path(seed, k): v for k, v in value.items() if k not in skip}
        details = [(child_path(seed, k), 'Supporting data: ' + k.replace('_', ' '))
                   for k in skip if k in value]
        if seed == '/research/maternal_lineage' and isinstance(value.get('results'), list):
            values.pop(seed + '/results', None)
            for i, candidate in enumerate(value['results']):
                base = seed + '/results/' + str(i)
                if not isinstance(candidate, dict):
                    values[base] = candidate
                    continue
                for key, item in candidate.items():
                    path = child_path(base, key)
                    if key == 'Range':
                        details.append((path, f'Candidate {i + 1}: assayed mitochondrial positions'))
                    else:
                        values[path] = item
        for path in ('/limits', '/conclusions/evidence_limits', '/conclusions/method'):
            from .bundle import at_path
            try:
                values[path] = at_path(report, path)
            except (KeyError, IndexError, TypeError):
                pass
        result.append({'seed': seed, 'label': label, 'area': area, 'details': details,
                       'records': values,
                       'selection_rule': 'Original result fields, status, method, sources and limitations. Only explicitly named bulk reference or coverage data has a separate optional route.'})

    panels = research.get('ancestry', [])
    if isinstance(panels, list):
        for i, panel in enumerate(panels):
            if isinstance(panel, dict):
                add('/research/ancestry/' + str(i), panel, 'ancestry',
                    str(panel.get('name') or f'Ancestry comparison {i + 1}') + ': result and limitations', ('points',))
    for key, (area, label) in TITLES.items():
        value = research.get(key)
        if isinstance(value, dict):
            skip = ('rows',) if key == 'frequencies' else ()
            add('/research/' + key, value, area, label, skip)
    return result


def path_label(path):
    """Describe storage locations without exposing raw JSON syntax as the label."""
    parts = path.strip('/').split('/')
    names = {'research': '', 'ancestry': 'Ancestry comparison', 'points': 'plotted reference sample',
             'paternal_lineage': 'Paternal lineage', 'maternal_lineage': 'Maternal lineage',
             'results': 'candidate', 'candidates': 'candidate', 'Range': 'assayed positions',
             'frequencies': 'Reference frequencies', 'rows': 'record', 'path': 'supporting branch',
             'score_coverage': 'Measured score coverage', 'imputed_score': 'Estimated score',
             'phase': 'Targeted phase', 'array_quality': 'Array quality'}
    labels = []
    for part in parts:
        if part.isdigit():
            labels.append(str(int(part) + 1))
        elif part in names:
            if names[part]:
                labels.append(names[part])
        else:
            labels.append(re.sub(r'[_~]', ' ', part))
    return ' · '.join(labels) or 'Report structure'
