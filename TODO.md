# Shared `nf-core/modules` follow-up

The pipeline standardizes gene identifiers and `gene_symbol` annotations in its AnnData outputs. When updating the shared `anndata/convert` module, add regression assertions that the conversion preserves AnnData `var_names` as Seurat feature names and SingleCellExperiment row names, and preserves `var["gene_symbol"]` as Seurat `meta.features$gene_symbol` and SingleCellExperiment `rowData()$gene_symbol`. Make those shared-module changes in `nf-core/modules` rather than editing the vendored module copy in this pipeline.

# test pipeline

# update nf-tests
