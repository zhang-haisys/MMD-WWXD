# External helper functions and references

The submitted scripts call helper files not included in the supplied directory. Examples found during path review:

- `/opt/scripts_repo/bin/`: `pheatmap1.R`, `plot_boxplot.R`, `plot_vioplot.R`, `差异分析.R`, `差异分析中文报告.R`, `plot_Vennplot.R`, enrichment plotting and enrichment functions, `BestCC.r`, survival helpers, and other utilities.
- `/work/run/wangtc/predata/ML_FeatureSelect.r`
- `/JJGene/data/xuwj/prepare/Gene_expression.R`
- External reference files under `/work/public/ceRNA/` and `/work/run/leiyp/...`.

These helper implementations and their licenses/provenance were not provided. Scripts that call them will not run as-is until the exact dependencies are obtained and configured. This repository does not fabricate replacement helper functions.
