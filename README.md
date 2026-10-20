# DNA Research Report

Turn a raw DNA text file into a research report using public genetic references. The app runs in your terminal, checks its findings and saves a report you can read or give to an LLM for a plain-language explanation.

Inspired by Promethease’s raw-DNA-to-report workflow, it generates its own findings without needing a Promethease report or an LLM to run.

For example, when the necessary markers are available, it can identify APOE patterns associated with late-onset Alzheimer’s risk. APOE ε2 is associated with lower risk in some populations, but does not guarantee protection. See the [National Institute on Aging’s explanation](https://www.nia.nih.gov/research/milestones/disease-mechanisms/disease-mechanisms-apoe-physiology).

## Open the app

Open the project folder in a terminal and run:

```sh
./start
```

**First launch takes longer.** Setup installs Python, analysis packages, Java and PLINK, then downloads public genetic references and builds the local databases. Allow at least 10 GiB of free space and tens of minutes. Setup does not read your DNA. If interrupted, run `./start` again; completed downloads are kept.

**Later launches reuse that setup.** The menus open quickly and analysis can run offline. Generating a report still takes processing time, shown by a progress bar.

Automatic setup supports macOS on Apple Silicon or Intel and Linux x86_64 with glibc. The full workflow has been tested on macOS. See [setup and recovery](docs/REFERENCE_SETUP.md) for details.

```text
DNA Research Report

What would you like to do?
  > Analyze a DNA file
    Export a saved report
    Try an invented example
    Quit
```

Use **Up/Down** and **Enter** to select. Choose **Back** to revisit a step, and scroll up to read earlier screens. You can browse to files and folders or paste a path; no analysis commands are needed.

1. Choose **Analyze a DNA file** and select your file.
2. Choose where to save the report. The app detects the file’s reference genome and DNA letter format, asking for help only if needed.
3. Review the settings and select **Generate report**.
4. The report folder is saved automatically. Select **Finish**, then open `report-bundle/overview.md` in your output folder. Other export formats are optional.

Reports default to a dated folder under `~/dna-reports/`. They must stay outside the project, and existing files are never overwritten. If the app cannot detect a setting, check your provider’s file instructions rather than guessing. Limited analysis requires your explicit choice.

**Export a saved report** creates another export without rerunning the analysis. **Try an invented example** lets you explore the flow without personal DNA. For exports or a limited three-rule example without full setup, use `./start --app-only`.

## Explain your report with an LLM

The flow is **overview and starting map → area or gene range → selected subject → required reading files**. Each reading plan names the files needed for that selection, including supporting evidence and warnings. Medication summaries also include other findings and unresolved calls assigned to the same gene. Research areas such as ancestry and lineages open with their results, methods and limitations; bulk plotting data and repeated coverage lists have separate optional routes.

Upload only the requested branch, not the entire folder or compressed archive. Files accumulate in a conversation, so exploring many subjects can still fill an LLM’s context.

1. Upload **`overview.md` and `index.md`** from `report-bundle/` to your chosen LLM.
2. Describe what you want to understand, or ask it to help you choose where to start.
3. Upload the requested maps until you reach a subject and its required reading plan.
4. Upload every file in that plan, then stop. Other subjects are optional; do not follow every link in the report.

Send complete files so their definitions, sources and warnings stay attached. A link to a file on your computer does not upload it. An assistant with access to the report folder can find the files itself.

You can start with this prompt:

> Help me explore this DNA report. Use the overview and index to suggest where to start, then ask for the files you need. Choose the relevant subject, request every file in its required reading plan, then stop. Do not expand into other subjects unless I ask. Explain only the findings in files I provide, in plain language, and cite their sources. Keep warnings, conflicting evidence and missing information clear. Distinguish DNA read from my file from estimates. Do not treat an association as a diagnosis or assume you have read the rest of the report.

The app does not connect to an LLM. Uploading files to an online service shares the genetic information in those files with that service.

## What the report includes

The export keeps all results produced by the analysis, including:

- Findings based on individual DNA markers or combinations of markers, with explanations and sources.
- Clinical database annotations, research study links, conflicting evidence and warnings.
- Ancestry comparisons, DNA-pattern frequencies, and maternal and paternal lineage estimates.
- Research scores, including an estimated BMI score, and a targeted analysis of how TPMT markers may be inherited together.
- What could be checked, what remains unresolved, and which results use estimates or incomplete data.

The output folder contains:

| File or folder | Purpose |
| --- | --- |
| `overview.md` | Start here for coverage and a guide to the report. |
| `index.md` | Choose an area to explore. |
| `maps/` | Area and gene-range maps, subject selectors and required reading plans. |
| `focus/` | Selected findings and research results with their supporting records and limitations. |
| `groups/` | Complete source records, including results outside the focused summaries. |
| `index.tsv` | Searchable list for finding specific genes or markers. |
| `full-report.json.gz` | Compressed copy of the complete report for storage and future exports. |
| `manifest.json` | File inventory and integrity checks. |

The same finding can be reached through several maps. Links use explicit marker IDs and source relationships; sharing a marker does not establish the same allele or medical meaning. Map labels help you navigate, not interpret a finding. Unclassified data remains available.

An assistant with local file tools can also decompress and search `full-report.json.gz` for a specific gene or marker, then read its supporting records and warnings. It does not need to load the entire report into a conversation.

To create another export, select `full-report.json.gz` through **Export a saved report**. Single-file JSON, Markdown and text exports are also available. See [export details and optional commands](docs/CLI.md).

## How it works

Python runs the analysis using local SQLite databases. Findings come from explicit rules, cited source text and report templates.

| Analysis | Tools and public references |
| --- | --- |
| Marker explanations and combinations | SNPedia via the public SNP Browser snapshot and SNPedia bot API |
| DNA positions and matching | Ensembl GRCh37, UCSC assembly chains, ClinGen Allele Registry and publications |
| Clinical evidence and test warnings | ClinVar and GenomePrep |
| Research study links | GWAS Catalog, with up to 2,000 study overlaps retained separately from personal effect estimates |
| Ancestry and frequencies | PLINK, NumPy and scikit-learn with 1000 Genomes and HGDP panels |
| Maternal and paternal lineages | HaploGrep 3 with Phylotree, and FamilyTreeDNA’s public Y tree |
| Scores and estimated DNA patterns | PGS Catalog and Beagle 5.5 with public reference haplotypes |

The app checks reference files for changes before analysis. [Methods](docs/METHODS.md) explains the calculations, matching rules and validation limits.

## Scope and privacy

This is a research tool, not a diagnosis, medication guide or complete assessment of genetic health. Missing or unsupported DNA markers can leave findings unresolved. Ancestry comparisons are not ancestry percentages, and research scores are not personal disease probabilities.

It does not reproduce every Promethease finding. Comprehensive drug-response analysis, dedicated blood-group reports, comparisons between people and a complete index by disease are outside its supported scope.

The project contains code, methods and independently invented examples. Keep personal DNA, reports and reference caches outside the repository. SNPedia content retains attribution under [CC BY-NC-SA 3.0 US](https://creativecommons.org/licenses/by-nc-sa/3.0/us/). See [source terms](docs/REFERENCE_SETUP.md#source-terms) before redistributing data or derived text.

## Development

```sh
python -m unittest discover -s tests -v
python tools/check_release.py
python -m tools.validate_indels --out /path/to/normalization-check.json
```

Tests cover DNA parsing and matching, missing data, warnings, reference integrity, terminal navigation, progress, cancellation and complete exports. Export tests reconstruct the original data from the report sections and check bounded maps, finite reading plans and warning selection. Focused extracts must retain the exact source fields and exclude unrelated records that happen to share a file or a marker with another rule. These checks do not establish an external LLM’s comprehension or medical accuracy. The normalization tool checks synthetic DNA insertions and deletions against bcftools. Earlier experiments are in `research/`.
