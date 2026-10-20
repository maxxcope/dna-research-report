import argparse
import json
from pathlib import Path
import sqlite3
import sys
from .engine import analyze
from .references import import_clinvar, import_gwas, fetch_clinvar
from .report import export_report, read_report
from .local import analyze_local, DEFAULT_REFS, read_input
from .progress import Progress
from .terminal import guided, QuitFlow, Terminal, describe_result, application
from .extras import cancel_commands

def write_report(report, out, formats=('bundle',)):
    formats = tuple(formats)
    if 'bundle' not in formats:
        formats = tuple(dict.fromkeys(('json',) + formats))
    paths = export_report(report, out, formats)
    for path in paths.values():
        print('Report written:', path)
    return paths


def selected_formats(value):
    return ('bundle', 'json', 'markdown', 'text') if value == 'all' else (value,)


def main(argv=None):
    parser=argparse.ArgumentParser(description='Local raw-DNA research reports. No Promethease required.')
    parser.add_argument('--exports-only',action='store_true',help=argparse.SUPPRESS)
    commands=parser.add_subparsers(dest='command')
    run = commands.add_parser('run', help='Guided terminal analysis with automatic build detection')
    run.add_argument('input', type=Path, nargs='?')
    run.add_argument('--refs', type=Path, default=DEFAULT_REFS)
    run.add_argument('--out', type=Path)
    a=commands.add_parser('analyze',help='Create a report from a raw text file')
    a.add_argument('input',type=Path);a.add_argument('--build',choices=['auto','37','38'],default='auto');a.add_argument('--refs',type=Path,default=DEFAULT_REFS)
    a.add_argument('--forward-strand',action='store_true',required=True,help='Confirm forward orientation from provider metadata')
    a.add_argument('--format', choices=['bundle','json','markdown','text','all'], default='bundle')
    a.add_argument('--progress', choices=['auto','plain','json'], default='auto', help='Real stage progress on stderr; JSON events are suitable for scripts')
    a.add_argument('--traits-only',action='store_true');a.add_argument('--out',type=Path,required=True)
    a.add_argument('--profile',type=Path,help='Full pipeline profile (default: REFS/pipeline.json when present)')
    a.add_argument('--work',type=Path,help='Empty directory to retain intermediate research outputs')
    a.add_argument('--basic',action='store_true',help='Explicitly run only the earlier ClinVar/traits report')
    k=commands.add_parser('import-snpedia',help='Prepare public SNPedia knowledge and rules')
    k.add_argument('--database',type=Path,required=True);k.add_argument('--criteria',type=Path,required=True)
    k.add_argument('--chain',type=Path);k.add_argument('--refs',type=Path,default=DEFAULT_REFS)
    k.add_argument('--coordinate-reference',type=Path,help='Optional prepared ClinVar catalogue for unique SNV positions')
    k.add_argument('--placement-reference',type=Path,help='Optional complete public Ensembl placement cache')
    demo=commands.add_parser('demo',help='Generate an entirely synthetic report');demo.add_argument('--out',type=Path,default=Path('../dna-reports/demo'));demo.add_argument('--refs',type=Path,default=DEFAULT_REFS)
    demo.add_argument('--format', choices=['bundle','json','markdown','text','all'], default='bundle')
    export = commands.add_parser('export', help='Export an existing complete JSON report without rerunning analysis')
    export.add_argument('report', type=Path);export.add_argument('--out', type=Path, required=True)
    export.add_argument('--format', choices=['bundle','json','markdown','text','all'], default='bundle')
    c=commands.add_parser('import-clinvar',help='Prepare a local public ClinVar VCF')
    c.add_argument('source',type=Path);c.add_argument('--build',choices=['37','38'],required=True);c.add_argument('--refs',type=Path,default=DEFAULT_REFS)
    g=commands.add_parser('import-gwas',help='Prepare a downloaded public GWAS Catalog TSV or ZIP')
    g.add_argument('source',type=Path);g.add_argument('--refs',type=Path,default=DEFAULT_REFS)
    f=commands.add_parser('fetch-clinvar',help='Download a public reference, without reading DNA')
    f.add_argument('--build',choices=['37','38'],required=True);f.add_argument('--refs',type=Path,default=DEFAULT_REFS)
    args=parser.parse_args(argv)
    try:
        if args.command in (None, 'run'):
            if args.command is None and not sys.stdin.isatty():
                parser.print_help();return 2
            if args.command is None:
                return application(limited=args.exports_only)
            return guided(getattr(args, 'refs', DEFAULT_REFS), getattr(args, 'out', None), getattr(args, 'input', None)) or 0
        elif args.command == 'export':
            report = read_report(args.report)
            for path in export_report(report, args.out, selected_formats(args.format)).values():
                print('Report written:', path)
        elif args.command=='analyze':
            data = read_input(args.input)
            profile = args.profile or args.refs / 'pipeline.json'
            with Progress(mode=args.progress) as progress:
                if not args.traits_only and not args.basic:
                    from .pipeline import full_report
                    if not profile.exists():
                        raise ValueError('Full pipeline profile missing. Prepare references, supply --profile, or explicitly choose --basic.')
                    report = full_report(data, profile, args.build, work=args.work, events=progress.event)
                else:
                    progress.event(dict(kind='detail', message='Running explicitly selected limited analysis'))
                    report = analyze_local(data, args.build, args.refs, args.traits_only)
            describe_result(report, Terminal())
            write_report(report, args.out, selected_formats(args.format))
        elif args.command=='import-snpedia':
            from .knowledge import import_knowledge
            print(json.dumps(import_knowledge(args.database,args.criteria,args.chain,args.refs/'snpedia.sqlite',args.coordinate_reference,args.placement_reference),indent=2))
        elif args.command=='demo':
            from importlib.resources import files
            data=files('genetics_app').joinpath('data/demo_dna.txt').read_bytes()
            if (args.refs/'pipeline.json').exists():
                from .pipeline import full_report
                with Progress() as progress:
                    report=full_report(data,args.refs/'pipeline.json',events=progress.event)
            else:report=analyze(data,traits_only=True)
            describe_result(report, Terminal())
            write_report(report,args.out,selected_formats(args.format))
        elif args.command=='import-clinvar':
            m=import_clinvar(args.source,args.refs/f'clinvar{args.build}.sqlite',args.build)
            print('Imported public ClinVar records:',m['records'])
        elif args.command=='import-gwas':
            m=import_gwas(args.source,args.refs/'gwas.sqlite');print('Imported public GWAS records:',m['records'])
        elif args.command=='fetch-clinvar':
            m=fetch_clinvar(args.build,args.refs/f'clinvar{args.build}.sqlite');print('Imported public ClinVar records:',m['records'])
    except QuitFlow:
        return 0
    except KeyboardInterrupt:
        cancel_commands()
        print('\nCancelled. Temporary calculation files have been removed.', file=sys.stderr);return 130
    except (OSError,ValueError,sqlite3.Error) as e:
        print('Error:',str(e),file=sys.stderr);return 2
    return 0

if __name__=='__main__':raise SystemExit(main())
