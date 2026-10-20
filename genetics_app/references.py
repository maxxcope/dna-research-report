"""Public-only reference preparation. Analysis never downloads data."""
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sqlite3
import tempfile
from datetime import datetime, timezone
import urllib.request
import zipfile

SCHEMA = 1

def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def connect_read(path):
    db=sqlite3.connect(Path(path).resolve().as_uri()+'?mode=ro',uri=True)
    db.row_factory=sqlite3.Row
    meta=dict(db.execute('select key,value from metadata'))
    if meta.get('schema')!=str(SCHEMA):
        db.close();raise ValueError('Unsupported reference database schema.')
    return db,meta

def metadata(path):
    db,m=connect_read(path);db.close();return m

def import_clinvar(source,destination,build,source_url='',progress=lambda _:None):
    source,destination=Path(source),Path(destination)
    if source.resolve()==destination.resolve():raise ValueError('Source and destination must differ.')
    destination.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent) as tmp:
        db=sqlite3.connect(str(Path(tmp)/'reference.sqlite'))
        db.execute('create table metadata(key text primary key,value text)')
        db.execute('create table variants(chrom text,pos integer,id text,ref text,alt text,rs text,gene text,classification text,review text,conditions text,origin text,conflicts text)')
        count=0;declared=None;date='unknown';batch=[]
        opener=gzip.open if source.suffix=='.gz' else open
        with opener(source,'rt',encoding='utf-8') as f:
            for line in f:
                if line.startswith('##reference='):
                    match=re.search(r'GRCh(37|38)',line)
                    declared=match.group(1) if match else None
                if line.startswith('##fileDate='):date=line.strip().split('=',1)[1]
                if line.startswith('#'):continue
                if declared != build:raise ValueError('ClinVar reference header does not match the requested build.')
                parts=line.rstrip('\n').split('\t')
                if len(parts)<8:raise ValueError('Malformed ClinVar VCF.')
                ch,pos,vid,ref,alt=parts[:5]
                info=dict(x.split('=',1) for x in parts[7].split(';') if '=' in x)
                ch=ch.removeprefix('chr')
                batch.append((ch,int(pos),vid,ref,alt,info.get('RS',''),info.get('GENEINFO',''),info.get('CLNSIG',''),info.get('CLNREVSTAT',''),info.get('CLNDN',''),info.get('ORIGIN',''),info.get('CLNSIGCONF','')))
                count+=1
                if len(batch)>=10000:
                    db.executemany('insert into variants values(?,?,?,?,?,?,?,?,?,?,?,?)',batch);batch=[]
                    if count % 500000 == 0:progress(f'ClinVar records indexed: {count:,}')
        if not count or declared!=build:raise ValueError('Empty or incompatible ClinVar reference.')
        db.executemany('insert into variants values(?,?,?,?,?,?,?,?,?,?,?,?)',batch)
        db.execute('create index locus on variants(chrom,pos)')
        meta=dict(schema=str(SCHEMA),kind='ClinVar',build=build,release_date=date,source_sha256=digest(source),source_url=source_url or 'local public ClinVar VCF import',records=str(count),prepared_utc=datetime.now(timezone.utc).isoformat())
        db.executemany('insert into metadata values(?,?)',meta.items());db.commit();db.close()
        Path(tmp,'reference.sqlite').replace(destination)
    progress(f'ClinVar index complete: {count:,} records')
    return meta

def import_gwas(source,destination,source_url='',progress=lambda _:None):
    source,destination=Path(source),Path(destination)
    if source.resolve()==destination.resolve():raise ValueError('Source and destination must differ.')
    destination.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent) as tmp:
        db=sqlite3.connect(str(Path(tmp)/'gwas.sqlite'))
        db.execute('create table metadata(key text primary key,value text)')
        db.execute('create table associations(rs text,trait text,pmid text,p real,effect text,effect_allele text,study text,population text,publication_date text,reported_gene text)')
        archive=None
        if source.suffix=='.zip':
            archive=zipfile.ZipFile(source)
            names=[n for n in archive.namelist() if not n.endswith('/') and n.endswith(('.tsv','.txt'))]
            if len(names)!=1:raise ValueError('Expected one TSV/TXT association file in GWAS ZIP.')
            stream=io.TextIOWrapper(archive.open(names[0]),encoding='utf-8-sig')
        else:stream=open(source,encoding='utf-8-sig')
        count=0;batch=[];latest=''
        with stream:
            reader=csv.DictReader(stream,delimiter='\t')
            if not {'SNPS','P-VALUE','PUBMEDID','DISEASE/TRAIT'}<=set(reader.fieldnames or []):raise ValueError('Unsupported GWAS Catalog table.')
            for r in reader:
                latest=max(latest,r.get('DATE ADDED TO CATALOG',''))
                try:p=float(r['P-VALUE'])
                except ValueError:continue
                if not math.isfinite(p) or not 0<p<=5e-8:continue
                ids=set(re.findall(r'\brs\d+\b',r['SNPS']+' '+r.get('STRONGEST SNP-RISK ALLELE','')))
                if r.get('SNP_ID_CURRENT','').isdigit():ids.add('rs'+r['SNP_ID_CURRENT'])
                for rs in ids:
                    batch.append((rs,r['DISEASE/TRAIT'],r['PUBMEDID'],p,r.get('OR or BETA',''),r.get('STRONGEST SNP-RISK ALLELE',''),r.get('STUDY ACCESSION',''),r.get('INITIAL SAMPLE SIZE',''),r.get('DATE',''),r.get('MAPPED_GENE','')))
                    count+=1
                if len(batch)>=10000:
                    db.executemany('insert into associations values(?,?,?,?,?,?,?,?,?,?)',batch);batch=[]
                    progress(f'GWAS annotations indexed: {count:,}')
        if archive:archive.close()
        if not count:raise ValueError('No eligible GWAS records found.')
        db.executemany('insert into associations values(?,?,?,?,?,?,?,?,?,?)',batch)
        db.execute('create index marker on associations(rs)')
        meta=dict(schema=str(SCHEMA),kind='GWAS Catalog',records=str(count),latest_entry_date=latest,source_sha256=digest(source),source_url=source_url or 'local public GWAS Catalog import',prepared_utc=datetime.now(timezone.utc).isoformat(),threshold='p <= 5e-8',matching='rsID literature overlap only; no effect-allele interpretation')
        db.executemany('insert into metadata values(?,?)',meta.items());db.commit();db.close()
        Path(tmp,'gwas.sqlite').replace(destination)
    progress(f'GWAS index complete: {count:,} annotations')
    return meta

def fetch_clinvar(build,destination):
    if build not in {'37','38'}:raise ValueError('Unsupported build')
    url=f'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh{build}/clinvar.vcf.gz'
    destination=Path(destination);destination.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent) as tmp:
        p=Path(tmp)/'clinvar.vcf.gz'
        with urllib.request.urlopen(url,timeout=60) as src, p.open('wb') as out:
            while chunk:=src.read(1024*1024):out.write(chunk)
        return import_clinvar(p,destination,build,url)
