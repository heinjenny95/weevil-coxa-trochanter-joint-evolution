# Submission consistency audit

`build_proxy_mapping.R` creates an explicit specimen-to-proxy-tip table from
`specimen_key.csv`. The output distinguishes the sampled genus from the 15
taxonomic proxy tips used for phylogenetic comparative analyses and reports the
number of specimens aggregated into each proxy tip.

```text
Rscript build_proxy_mapping.R specimen_key.csv taxonomic_proxy_mapping.csv
```

The mapping is descriptive: it does not alter the phylogeny or reassign taxa.
Its purpose is to make the unavoidable taxonomic approximation auditable.

`refresh_shape_geometry_tables.R` rebuilds the concise main-dataset
shape--geometry summary from the canonical regression table and validates the
expected sample size (63). `sync_supplementary_source_data.py` then
synchronizes Supplementary Tables 13--39 with the robust point-estimate
outputs, labels geometry-dependent files as `main_dataset`, removes private
workstation paths, normalizes all released tables to comma-delimited CSV with decimal
points, and regenerates the source-data checksum manifest.

```text
Rscript refresh_shape_geometry_tables.R <repository-root>
python sync_supplementary_source_data.py <repository-root> <robust-point-estimates-root>
```
