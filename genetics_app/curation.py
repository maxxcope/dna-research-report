"""Reviewed public sequence definitions, independent of personal inputs."""
import hashlib
from importlib.resources import files
import json


def definitions():
    data=files('genetics_app').joinpath('data/source_metadata.json').read_bytes()
    return json.loads(data)['definitions'],hashlib.sha256(data).hexdigest()


def apply_definition(marker,source_text,source_alleles,definition):
    if marker['id']!=definition['id'] or hashlib.sha256(source_text.encode()).hexdigest()!=definition['source_content_sha256']:return marker
    source=definition['source_change']
    if len(source)!=2 or set(source_alleles)!=set(source):return marker
    placements=definition['placements'];orientations={}
    for build,p in placements.items():
        if not (len(p['ref'])==len(p['alt'])==1 and set(p['ref']+p['alt'])<=set('ACGT') and p['ref']!=p['alt']):return marker
        if marker['chrom'] and marker['chrom']!=p['chrom']:return marker
        if marker['positions'].get(build) not in (None,p['pos']):return marker
        direct=source==[p['ref'],p['alt']]
        reverse=[a.translate(str.maketrans('ACGT','TGCA')) for a in source]==[p['ref'],p['alt']]
        if direct==reverse:return marker
        orientations[build]='plus' if direct else 'minus'
        if marker['orientation'] and marker['orientation']!=orientations[build]:return marker
    marker['curated_snv_definition']=definition
    marker['chrom']=placements['37']['chrom'];marker['gene']=marker['gene'] or definition['gene']
    marker['positions'].update({b:p['pos'] for b,p in placements.items()})
    marker['orientations_by_build']=orientations
    return marker
