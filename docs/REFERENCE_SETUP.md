# Reference setup

Public reference preparation is separate from DNA analysis. Setup has no genotype-file argument; analysis reads prepared files locally. Keep the reference cache and personal reports outside the repository.

## Automatic setup

Run `./start` from the project folder. No Python, Java, PLINK, package-manager installation or environment activation is needed first. The launcher downloads a version-pinned [uv](https://docs.astral.sh/uv/) binary, verifies its SHA-256, uses uv to manage Python 3.12.15, and installs the exact packages and hashes in `requirements.lock`. Java 21 from [Eclipse Temurin](https://adoptium.net/temurin/releases/) and [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/) are downloaded from their official releases and checked against the pinned inventory in `genetics_app/data/install_tools.json` before execution. The pinned macOS PLINK archive contains both Intel and Apple Silicon binaries.

Automatic installation targets macOS Apple Silicon/Intel and Linux x86_64 with glibc. First launch needs `curl`, `tar`, a SHA-256 utility and trusted HTTPS certificates, normally provided by these systems. Native Windows, Linux ARM and Alpine/musl are not supported by this launcher. Windows users can use an x86_64 Linux distribution under WSL2; that route and Linux have not been exercised end to end here. Operating-system security restrictions still apply. Setup never runs `sudo`, accepts OS licence prompts or modifies system Python or your shell configuration.

The default installation folder is `~/.local/share/dna-research-report/` (or `$XDG_DATA_HOME/dna-research-report/`). It holds managed Python, the analysis environment, tools, downloads, references and `setup.log`. Set `DNA_RESEARCH_HOME` to a different folder before starting if needed. Set `DNA_RESEARCH_REFS` to adopt an existing complete reference cache. Both must be outside the project. Reports remain separately under `~/dna-reports/` or your selected destination.

Allow at least **10 GiB free** for setup, reference indexes and calculation work. Downloads and index construction can take tens of minutes. The launcher checks free space before full reference preparation. Actual requirements depend on source releases and report sizes; disk exhaustion is reported rather than treated as successful setup. Beagle uses up to two workers with a 1.5 GB Java heap per worker; Python and the ancestry panels need additional memory.

### Recovery and repeat launches

- **Interrupted setup:** run `./start` again. Completed files, rule/API batches and reference windows are reused. Partial HTTP downloads resume when the server provides a matching ETag or Last-Modified validator; otherwise that file restarts. An unfinished imputation window is rebuilt. Indexes and the final profile are published only when complete.
- **Connection or provider error:** bounded download retries run automatically. If the provider remains unavailable, the terminal offers Retry, opening the app for exports/examples, or Quit. Setup never labels an incomplete reference profile as ready.
- **Offline use:** after a successful setup, launches reuse the environment and references without contacting providers. Analysis reads local references. Exports need neither the original DNA nor another analysis. `./start --app-only` can open the menus before full references are ready, once the managed Python runtime is installed. This fallback needs no research packages and offers only exports and the three-rule example.
- **Damaged files:** cached checksums are checked, not silently replaced with new hashes. Move the named damaged file aside and retry, or select a new reference directory. Keep the old cache if it is needed to reproduce a prior report.
- **Two launches:** a process lock prevents concurrent installers from changing the same environment. A second lock protects a shared reference cache even when installers use different application folders. The OS releases it after a crash; there is no stale lock directory to delete. Analysis begins after this lock is released.
- **Moved folder or changed dependencies:** the launcher checks package versions/imports and repairs its disposable environment when needed. Moving a cache may require reference preparation again to update executable paths. Personal reports are not removed.

`./start --setup-only` prepares and verifies the installation without opening the menus, returning a nonzero exit code if setup fails. Subsequent runs still verify profile identity and file sizes/timestamps. Changed files trigger full hash verification, and every analysis independently verifies all profile hashes. This avoids reading gigabytes on every menu launch while preserving analysis integrity.

For a new public reference release, choose a new cache directory. The launcher does not automatically update a working reference snapshot. Download accessibility does not override [source terms](#source-terms).

## Manual installation

This is an alternative for an existing Python environment or a system that needs custom tools:

```sh
python -m pip install '.[research]'
python -m genetics_app.setup --refs /path/to/genetics-references \
  --plink /path/to/plink --java /path/to/java
python -m genetics_app run --refs /path/to/genetics-references
```

Use PLINK **1.9**, not PLINK 2, and Java 8 or newer. `DNA_RESEARCH_PLINK` and `DNA_RESEARCH_JAVA` can override managed tool downloads with existing working executables. Executables must match the OS/architecture; current official Linux PLINK builds require glibc 2.31 or newer. A compatible custom executable can be used on an older system. Tool overrides are tested before reference preparation.

## Complete GRCh37 profile

Setup prepares:

| Reference | Used for |
| --- | --- |
| [Public SNP Browser snapshot](https://github.com/jonluca/snp-browser), plus [SNPedia’s bot API](https://bots.snpedia.com/api.php) | Genotype text, source metadata and public genoset rules |
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/) | Independent variant classifications |
| [Ensembl GRCh37 variation API](https://grch37.rest.ensembl.org/documentation/info/variation_post) | Public SNV coordinates, allele definitions and documented rsID aliases |
| [GWAS Catalog](https://www.ebi.ac.uk/gwas/docs/file-downloads/) | Literature leads |
| [GenomePrep](https://zenodo.org/records/5047472) | Array-content clusters and quality flags |
| [TCAG 1000 Genomes](https://tcag.ca/tools/1000genomes.html), [HGDP](https://zenodo.org/records/17375732) | Ancestry comparison and reference frequencies |
| [FamilyTreeDNA public Y tree](https://www.familytreedna.com/public/y-dna-haplotree/) | Paternal branch comparison |
| [HaploGrep 3](https://haplogrep.readthedocs.io/en/latest/installation/) and Phylotree rCRS 17.3 | Maternal lineage |
| [PGS Catalog](https://www.pgscatalog.org/downloads/) | BMI, CAD and height score definitions |
| [Beagle 5.5 and public 1000 Genomes reference haplotypes](https://faculty.washington.edu/browning/beagle/beagle.html) | Local imputation and targeted phase |
| [UCSC hg19-to-hg38 chain](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/) | Public coordinate reconciliation |

Reference windows cover every target in the public BMI score, plus a fixed TPMT region. Their selection does not depend on a person's DNA. The imputation panel uses 503 European reference samples; applicability outside that ancestry is limited.

Ensembl requests are selected from gaps in the entire public source catalogue and the public rule inventory. Setup checkpoints batches of at most 200 IDs with three workers. It records the source, criteria and chain hashes and rejects incomplete or incompatible caches. No personal marker list is used for these requests; DNA analysis makes no Ensembl requests.

The cache contains `pipeline.json`, the analysis profile, and `downloads.json`, the URL/hash inventory. Setup resumes cached downloads and rejects known checksum changes. Prepared source and rule hashes are recorded; runtime verifies the files listed in the profile. Data paths are relative to the cache, while installed executable paths remain explicit. For a new reference release, prepare a new directory rather than silently mixing old and new data.

## Run with the cache

```sh
python -m genetics_app run --refs /path/to/genetics-references
python -m genetics_app analyze /path/to/dna.txt --forward-strand \
  --refs /path/to/genetics-references --out /path/to/private-report
```

The main terminal menu provides file/folder browsing, automatic build detection, matching reference selection and Back navigation. Choosing **Export a saved report** opens existing JSON through the same menus. CLI `--work /path/to/empty-directory` retains intermediate results; without it, temporary files are removed after the run. A retained work directory must be empty to prevent reuse of another input’s outputs.

## Input formats and GRCh38

Supported inputs are uncompressed UTF-8 four-column marker/chromosome/position/genotype files, including CSV, and five-column files with separate alleles. The parser accepts GRCh37 and GRCh38 forward-strand exports. It rejects conflicting build declarations, reverse-strand exports, VCF inputs and ambiguous PAR labels. Ambiguous indel codes remain unresolved; supported literal codes and simple heterozygous sequence mappings follow the guards in [Methods](METHODS.md).

The terminal app reads explicit forward/plus descriptions in the file's comment header and proceeds without a separate confirmation when the format is established. If none is found, **I'm not sure** explains where to check; confirmation is required before proceeding. A provider name or build number alone does not establish this format. Explicit reverse/minus or conflicting direction labels are rejected. For background, see [Ancestry's description of raw-data letters and strand](https://help.ancestry.com/hc/en-us/articles/53933317283603-Downloading-DNA-Data).

The complete default research profile uses GRCh37. A GRCh38 file needs matching references; the GRCh37 profile cannot be used by changing the header. GRCh38 currently supports the basic ClinVar/curated-trait workflow:

```sh
python -m genetics_app fetch-clinvar --build 38 --refs /path/to/grch38-references
python -m genetics_app analyze /path/to/grch38-dna.txt --forward-strand --basic \
  --refs /path/to/grch38-references --out /path/to/private-report
```

Existing public ClinVar VCF and GWAS TSV/ZIP files can be imported with `import-clinvar` and `import-gwas`. The `import-snpedia` command accepts a public snapshot, complete criteria inventory and optional chain. Use `python -m genetics_app COMMAND --help` for arguments.

## Source terms

The repository contains application code and synthetic fixtures; reference corpora are downloaded separately. SNPedia content remains attributed to its contributors under CC BY-NC-SA 3.0 US, including in generated reports. Check each provider's terms before redistributing reference data or derived source text. An accessible download is not a blanket commercial-use license.
