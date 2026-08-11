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
