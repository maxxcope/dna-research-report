# Advanced commands

The normal interface is the terminal menu opened by `python -m genetics_app`. These commands are optional for scripts, automation or repeatable configurations.

## Analyze without prompts

```sh
python -m genetics_app analyze /path/to/dna.txt --forward-strand \
  --refs /path/to/references --out /path/to/private-report
```

The default `--format bundle` writes `report-bundle/`, containing an overview, layered navigation maps, a complete TSV lookup, reading files and `full-report.json.gz`. No uncompressed duplicate is created by default.

Choose `json`, `markdown` or `text` for a legacy single-file export; analysis also saves JSON with a selected legacy readable format. `all` writes the bundle and all three single-file formats and can use substantially more disk space. Progress goes to stderr; `--progress json` emits actual stage and work-count events. `--build auto` is the default. A manual build cannot contradict a declared header.

`--basic` explicitly selects the limited ClinVar/three-trait workflow. `--traits-only` selects only the small curated panel. Neither replaces the full reference profile.

## Export an existing report

```sh
python -m genetics_app export /path/to/report-bundle/full-report.json.gz \
  --out /path/to/another-folder
```

No DNA file or reference cache is needed. Existing JSON reports and bundle directory paths are also accepted by this command. Use `--format json`, `markdown`, `text` or `all` when needed. Existing output files and bundle directories, including empty directories, are never overwritten. The same action is available through **Export a saved report** in the main menu; select the JSON file or gzip archive.

For an LLM conversation, supply `overview.md` and `index.md`. The LLM requests an area map, follows smaller subject/range maps when needed, then requests the reading files. Choose a subject and follow its required reading plan, including all its sections and warnings; then stop. A subject selector is a menu, not a checklist. Other subjects are optional. Upload the complete files, including definitions and qualifications. Map labels are source previews, not interpretations; local links do not give a chat access to the files.

`index.tsv` remains a complete path-to-file lookup for local searching and programmatic use. It can be too large to upload and is not needed for the normal map workflow.

All report destinations must be outside the application folder so personal results cannot be included in the public project by accident.

## Report file structure

```text
report-bundle/
  overview.md
  index.md
  maps/
    associations/index.md
    clinical/index.md
    ...
    subjects/
    companions/
    focus/
  focus/
    <subject>-0001.md
    ...
  groups/
    findings/part-0001.md
    ...
  index.tsv
  full-report.json.gz
  manifest.json
```

Maps have at most 24 links and 6,000 UTF-8 bytes. Large maps split into smaller ranges automatically; no exhaustive part list is placed in the starting map. Available areas are generated from the report, and findings can appear on multiple routes. Unknown fields remain in **Other retained data**. These routes are not a comprehensive disease taxonomy.

Existing evidence summaries and research results have focused reading sets in `focus/`. These retain the selected summary, its forward marker/rule references, saved clinical context and linked warning records. Gene-labelled medication summaries also include other directly annotated records for that gene, including unresolved calls; a matched pattern must not hide another measured allele with a different source label. They do not pull in other rules merely because those rules use the same marker. A plan lists all required sections and an explicit stopping point. It is a selection of saved evidence, not a complete clinical phenotype or all biological knowledge about a gene.

Ancestry result plans keep the original population comparisons, marker coverage, validation and limitations together. Plotted reference points have a separate optional route. Maternal lineage plans retain every candidate and its quality and mismatch fields; repeated lists of assayed positions have optional routes. Paternal candidates retain their supporting paths and conflicts. Reference frequency overviews retain their scope and limitations, with detailed rows available separately. Unknown result fields stay in the result plan by default. These are selections of original fields, not new scientific interpretations.

Other source records remain in `groups/`, packed by explicit primary marker ID. Their subject plans name the exact record paths to read; unrelated records sharing a file are outside the selection. Companion maps are subject selectors, not instructions to explore every subject. Explicit probe aliases connect warnings. A link does not assert identical alleles or clinical meaning, and source prose is not used to invent relationships. Unmatched rules remain available in the complete source catalogue.

Reading files have at most 24,000 UTF-8 bytes by default. Each includes its own definitions and qualifications. Repeated text is defined once within a file, and tables name their columns once. Oversized records are split and linked; readers must collect all sections before interpreting them. Together, the files preserve every retained field. Byte budgets do not guarantee a fit in every model’s context.

`full-report.json.gz` preserves the complete structured analysis. Compression reduces storage size, not the amount of information an LLM must read after decompression. `manifest.json` records archive, reading-file, lookup and map checksums for integrity verification, plus the original record paths selected for each focused view. Focused views duplicate selected evidence for convenient reading; the canonical parts and archive still preserve every original field. Version 4 adds readable research-result plans and optional bulk-data routes; existing saved JSON and gzip reports can be re-exported without DNA or reference downloads.

## Search the complete report locally

The gzip archive remains a full report. An assistant with file tools can decompress and search it without sending every record to the model. Search results identify candidate subjects; choosing a subject determines which supporting records and warnings to read.

For Python tools, the same selection used by focused exports is available directly:

```python
from genetics_app.report import read_report
from genetics_app.retrieval import EvidenceIndex

index = EvidenceIndex(read_report("/path/to/full-report.json.gz"))
matches = index.search("CYP2C19")
for match in matches:
    print(match["path"], match["label"])

# Choose a relevant path from those results, rather than selecting every match.
chosen_path = input("Record path: ").strip()
selection = index.select([chosen_path])
print(selection["records"])
print(selection["stop"])
```

Search uses exact words in saved gene, marker and source-label metadata. It is not a disease-question interpreter. Selection follows explicit forward dependencies and warning aliases, not neighboring records or reverse links into other rules. Study associations remain separately searchable and are included when explicitly selected. `no_separate_record_for` means no standalone saved record was found for a reference; it does not mean a genotype or condition is absent.

## Set starting locations

```sh
python -m genetics_app run --refs /path/to/references --out /path/to/private-report
```

This still opens the guided terminal flow. `DNA_RESEARCH_REFS` can set the default cache location. Otherwise it is `~/.cache/dna-research-report/references`.

## Synthetic example

```sh
python -m genetics_app demo --out ../dna-demo
```

Without prepared references, this demonstrates three curated traits. With a complete profile, it runs the full pipeline and shows missing coverage. The menu’s **Try an invented example** provides the same behavior.
