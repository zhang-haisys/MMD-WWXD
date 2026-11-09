# Reproducibility and execution notes

This package organizes scripts that were supplied from an existing project. It does not assert that a fresh checkout can reproduce the paper without additional files or edits.

The scripts reference multiple historical Linux and Windows paths (for example `/JJGene/...`, `/work/...`, `/opt/scripts_repo/bin/...`, and `C:/Users/...`). They also assume intermediate results and directories already exist. Categorized scripts are copies for navigation; originals are retained unchanged in `archive/original/`.

Several analyses call helper functions that are not part of the supplied script set. Their names and references are recorded in `MISSING_HELPERS.md`. Do not replace these helpers with guessed implementations if exact reproduction is required.

Before reporting results from a rerun, record the R version, package versions, data accessions and download dates, preprocessing decisions, random seeds, and any local changes. The source archive contains no synthetic `sessionInfo()` or claims of successful end-to-end execution.
