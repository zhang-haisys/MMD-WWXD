# Input-data inventory and expected placement

The supplied scripts expect a project-specific `0.data/` directory and several generated intermediate files. The data themselves are intentionally not included in this GitHub package.

Expected inputs mentioned in the scripts include TCGA-LUAD expression/clinical/mutation data; GEO series matrices and phenotype tables for GSE11969, GSE202859 and external validation; drug-target files; TIDE results; and precomputed ceRNA/TF/RBP reference tables. Exact filenames vary between script versions. Inspect each script's `read.csv`, `read.table`, and related calls before running it.

Do not upload credentials, restricted datasets, or identifiable records. Public dataset accessions and source URLs should be recorded in the eventual study-specific data manifest.
