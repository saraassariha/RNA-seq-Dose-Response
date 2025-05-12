# RNA-seq Dose Response

Compare gene expression between low- and high-dose groups using edgeR. The analysis filters low-count genes, applies TMM normalization and fits a quasi-likelihood model to identify differentially expressed genes.

## Quick start

Requires R and the edgeR, limma and statmod packages. Tested with R 4.5.2.

Run from the repository root:

```sh
Rscript install_packages.R
Rscript run.R
```

## Input

[`data/counts.tsv`](data/counts.tsv) contains integer read counts, with gene identifiers in the first column and 15 sample columns: `Low_1`–`Low_5`, followed by `High_1`–`High_10`.

To use another dataset, update the input path, group definitions and sample-name check in [`analysis.R`](analysis.R). Supply raw counts rather than normalized expression values.

## Output

Results are saved in `results/`:

- `differential_expression.csv`: genes with absolute log2 fold change > 1 and FDR < 0.05.
- `top_genes.csv`: the first ten significant genes, ranked by FDR.
- `plots.pdf`: the mean log2-CPM distribution before filtering.

Positive log2 fold changes indicate higher expression in the high-dose group. The included dataset yields **499 significant genes**. Additional sample metadata is needed for biological interpretation.

## Acknowledgments

University of Skövde for the example dataset and learning resources.
