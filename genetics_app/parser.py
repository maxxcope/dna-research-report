"""Strict input parsing and conservative locus resolution."""
from dataclasses import dataclass
from collections import Counter, defaultdict
import csv
import hashlib
import io
import re

MAX_BYTES = 64 * 1024 * 1024
MAX_ROWS = 2_000_000

class InputError(ValueError):
    pass

@dataclass(frozen=True)
class Call:
    marker: str
    chrom: str
    pos: int
    gt: str

@dataclass
class Genome:
    build: str
    build_source: str
    sha256: str
    rows: list
    by_position: dict
    by_id: dict
    bad_ids: set
    qc: dict

    def locus(self, chrom, pos):
        rows = self.by_position.get((chrom, pos), [])
        if not rows:
            return 'not_assayed', '', []
        ids = sorted({r.marker for r in rows})
        if any(i in self.bad_ids for i in ids):
            return 'conflicting_marker_id', '', ids
        calls = {r.gt for r in rows if r.gt != '--'}
        if len(calls) > 1:
            return 'conflicting_probes', '', ids
        if not calls:
            return 'no_call', '', ids
        gt = next(iter(calls))
        if not set(gt) <= set('ACGT'):
            return 'indel_code_unresolved', gt, ids
        if chrom.isdigit() and len(gt) != 2:
            return 'unexpected_autosomal_ploidy', gt, ids
        return 'called', gt, ids

def chromosome(s):
    s = s.strip().upper().removeprefix('CHR')
    return {'23':'X', '24':'Y', '26':'MT', 'M':'MT'}.get(s, s)

def detect_build(data):
    try:
        text = data.decode('utf-8-sig')
    except UnicodeDecodeError as e:
        raise InputError('Provide an uncompressed UTF-8 DNA text/CSV file.') from e
    if '\x00' in text:
        raise InputError('Binary input is unsupported.')
    headers = '\n'.join(l for l in text.splitlines() if l.lstrip().startswith('#'))
    found = set(re.findall(r'(?:GRCh\s*|hg\s*|build\s*[:=]?\s*)(36|37|38|18|19)\b', headers, re.I))
    found = { {'19':'37','18':'36'}.get(x,x) for x in found }
    if len(found) > 1:
        raise InputError('Conflicting genome builds in the header. Resolve the source metadata first.')
    detected = next(iter(found), None)
    if detected and detected not in {'37','38'}:
        raise InputError('Only GRCh37 and GRCh38 are supported; this header names another build.')
    return detected


def detect_strand(data):
    """Read explicit export-direction metadata; never infer it from DNA calls."""
    try:
        text = data.decode('utf-8-sig')
    except UnicodeDecodeError as error:
        raise InputError('Provide an uncompressed UTF-8 DNA text/CSV file.') from error
    headers = [line.lstrip().lstrip('#').strip() for line in text.splitlines() if line.lstrip().startswith('#')]
    word = re.compile(r'\b(forward|plus|reverse|minus)\b(?:\s*\([+-]\))?(?:\s*-\s*|\s+)strand\b', re.I)
    field = re.compile(r'\b(?:strand|orientation)\s*[:=]\s*(forward|plus|reverse|minus|[+-])(?=\s|$|[.,;])', re.I)
    evidence = {'forward': [], 'reverse': []}
    for index, line in enumerate(headers):
        for pattern in (word, field):
            for match in pattern.finditer(line):
                previous = headers[index - 1] if index and re.match(r'^on\b', line, re.I) else ''
                context = previous + ' ' + line
                if pattern is word and not re.search(r'\b(?:genotypes?|alleles?|calls?|data|reported|oriented)\b', context, re.I) and not re.fullmatch(r'(?:forward|plus|reverse|minus)(?:\s*\([+-]\))?(?:\s*-\s*|\s+)strand\.?', line, re.I):
                    continue
                prefix = re.split(r'[,;.!?]', line[:match.start()])[-1]
                if re.search(r'\b(?:not|no|unknown|unconfirmed|unverified|whether|might|may|neither)\b', ' '.join(prefix.split()[-7:]), re.I):
                    continue
                if re.match(r'\s*(?:is\s+)?(?:not\s+(?:confirmed|verified|known)|unknown|unconfirmed|unverified)', line[match.end():], re.I):
                    continue
                direction = 'forward' if match[1].lower() in {'forward', 'plus', '+'} else 'reverse'
                if line not in evidence[direction]:
                    evidence[direction].append(line)
    if evidence['forward'] and evidence['reverse']:
        raise InputError('The file header names both forward/plus and reverse/minus data. Ask the provider for a consistently labelled export; do not guess.')
    direction = 'forward' if evidence['forward'] else 'reverse' if evidence['reverse'] else 'unknown'
    return direction, evidence.get(direction, [])


def parse_genome(data, build='auto', strand='forward'):
    if strand != 'forward':
        raise InputError('Only forward-strand exports are supported. Convert and verify the orientation first.')
    if len(data) > MAX_BYTES:
        raise InputError('Input exceeds the 64 MiB limit.')
    if detect_strand(data)[0] == 'reverse':
        raise InputError('This file labels its DNA letters as reverse/minus. Request a forward/plus export from the provider; changing the header does not convert the data.')
    detected = detect_build(data)
    text = data.decode('utf-8-sig')
    headers = '\n'.join(l for l in text.splitlines() if l.lstrip().startswith('#'))
    if build not in {'auto','37','38'}:
        raise InputError('Choose auto, GRCh37, or GRCh38.')
    if build != 'auto' and detected and detected != build:
        raise InputError('Selected build contradicts the input header.')
    resolved = detected if build == 'auto' else build
    if resolved is None:
        raise InputError('Genome build is absent. Select GRCh37 or GRCh38 from the provider metadata; do not guess.')
    rows = []
    fmt = None
    for line_no, line in enumerate(io.StringIO(text), 1):
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        fields = next(csv.reader([s])) if ',' in s else s.split()
        fields = [f.strip() for f in fields]
        if fields[0].lower() in {'rsid','rs_id','snp','marker','marker_id'} and len(fields)>2 and fields[2].lower() in {'position','pos'}:
            continue
        if len(fields) not in (4,5):
            raise InputError(f'Line {line_no}: expected four columns or five columns with two alleles. VCF is unsupported.')
        if fmt and fmt != len(fields):
            raise InputError(f'Line {line_no}: mixed four- and five-column formats.')
        fmt = len(fields)
        marker,ch,p = fields[:3]
        if not re.fullmatch(r'[A-Za-z0-9_.:-]{1,100}', marker):
            raise InputError(f'Line {line_no}: invalid marker identifier.')
        ch = chromosome(ch)
        if ch not in {str(i) for i in range(1,23)} | {'X','Y','MT'}:
            raise InputError(f'Line {line_no}: unsupported chromosome (including ambiguous PAR/25 labels).')
        try:
            pos = int(p)
        except ValueError as e:
            raise InputError(f'Line {line_no}: invalid position.') from e
        if not 1 <= pos <= 300_000_000:
            raise InputError(f'Line {line_no}: position is outside the accepted range.')
        gt = ''.join(fields[3:]).upper()
        if gt in {'--','-','00','0','NN','N','..','.'}:
            gt = '--'
        elif not re.fullmatch(r'[ACGTID]{1,2}', gt):
            raise InputError(f'Line {line_no}: unsupported or partial genotype call.')
        else:
            gt = ''.join(sorted(gt))
        rows.append(Call(marker,ch,pos,gt))
        if len(rows)>MAX_ROWS:
            raise InputError('Too many rows; this app accepts array exports, not whole-genome VCF.')
    if not rows:
        raise InputError('No genotype rows found.')
    bypos,byid = defaultdict(list),defaultdict(list)
    for r in rows:
        bypos[(r.chrom,r.pos)].append(r)
        byid[r.marker].append(r)
    bad = {i for i,rr in byid.items() if len({(r.chrom,r.pos,r.gt) for r in rr})>1}
    g = Genome(resolved,'header' if detected else 'user declaration',hashlib.sha256(data).hexdigest(),rows,dict(bypos),dict(byid),bad,{})
    states = Counter(g.locus(*k)[0] for k in bypos)
    g.qc = dict(rows=len(rows),unique_markers=len(byid),unique_positions=len(bypos),missing_rows=sum(r.gt=='--' for r in rows),
                provider='23andMe' if re.search(r'\b23andme\b',headers,re.I) else 'unspecified',
                duplicate_marker_ids=sum(len(rr)>1 for rr in byid.values()),conflicting_marker_ids=len(bad),
                duplicate_positions=sum(len(rr)>1 for rr in bypos.values()),locus_states=dict(states),
                synthetic_input_declared=bool(re.search(r'(?im)^#.*synthetic',headers)),
                chromosomes=dict(Counter(r.chrom for r in rows)),format='four-column' if fmt==4 else 'five-column')
    return g
