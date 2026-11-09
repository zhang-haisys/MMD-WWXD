# MMD + WXD LUAD analysis source code

This repository organizes the supplied R scripts for the MMD + WXD non-small-cell lung cancer study. It is a source-code archive and organization release; the scripts have not been validated as an end-to-end reproducible pipeline.

## Repository layout

- `scripts/01_preprocessing/`: TCGA/GEO retrieval, expression preprocessing, and batch adjustment.
- `scripts/02_differential_expression/`: resistance and tumor/control differential analyses.
- `scripts/03_gene_selection/`: intersections, enrichment, and drug-target processing.
- `scripts/04_modeling/`: machine learning, diagnostic models, subtype analysis, GSEA, and drug-sensitivity prediction.
- `scripts/05_followup/`: TMB, immune checkpoints, TIDE, interaction networks, pan-cancer, and external validation.
- `scripts/auxiliary/`: helper scripts supplied with the source package.
- `archive/original/`: unchanged copies of all supplied scripts, preserving original names and contents.
- `data/`: directory for user-provided inputs; data files are not included in this release.
- `config/`: local path configuration example.
- `docs/`: dependency, path, and reproducibility notes.

## Before running

1. Review `docs/REPRODUCIBILITY_NOTES.md` and `docs/INPUT_DATA.md`.
2. Obtain the required public datasets and place them in the documented local data layout.
3. Supply the custom helper functions and external reference files listed in `docs/MISSING_HELPERS.md`.
4. Configure local paths in `config/config.example.R` and adapt scripts as needed.

The supplied scripts contain machine-specific absolute paths, assumptions about pre-existing intermediate files, and custom functions from external script libraries. The `archive/original/` files are preserved for traceability. The categorized copies in `scripts/` are provided for navigation and may still require path/configuration edits before execution.

## Data and privacy

Do not commit controlled-access data, patient-level data, credentials, tokens, or large intermediate results. Public TCGA/GEO data should be retrieved under the respective repository terms and cited appropriately. This repository does not redistribute the source datasets.

## License

See `LICENSE` and `docs/THIRD_PARTY_NOTICES.md`. External helper code is not redistributed unless its license and provenance permit it.
