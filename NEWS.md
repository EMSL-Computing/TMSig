# TMSig 1.1.1 (2024-07-02)

**BUG FIXES**

- `cameraPR.matrix`: `data.table::setcolorder` arguments `before` and `after` were added in data.table v1.15.0 (2024-01-30). Rather than requiring a recent minimum version of data.table, the code has been modified to no longer depend on these arguments.


# TMSig 1.1.0 (2024-06-26)

**NEW FEATURES**

- `cameraPR.matrix`: added `alternative` parameter to allow users to perform one-sided tests. A warning will be issued if attempting to use a one-sided test with the parametric version of CAMERA-PR, as this is generally not recommended.
- `cameraPR.matrix`: added `min.size` parameter to allow users to specify a minimum set size for testing.
- `similarity` and `cluster_sets`: added support for Ōtsuka set similarity.
- Documentation for most functions has been updated.

**BUG FIXES**

- `cameraPR.matrix`: fixed variance calculation when `use.ranks=TRUE` and `inter.gene.cor=0`.
- `enrichmap`: fixed checks for duplicate contrast-set pairs. 
- `enrichmap`: if column or row labels are specified by the user through the `heatmap_args` parameter, the maximum text width or height is now determined from these new labels to avoid overlapping plot elements.
- `enrichmap`: default of `statistic_column` parameter changed to "TwoSampleT" to match output of `cameraPR.matrix`.


# TMSig 1.0.0 (2024-05-06)

**BREAKING CHANGES**

- Rename package from ostRich to TMSig: Tools for Molecular Signatures.
- Rename `prepare_gene_sets` vignette to `GO_analysis`. Add example code for `cameraPR.matrix`.

**NEW FEATURES**

- Add fast `cameraPR.matrix` method based on `limma::cameraPR.default`. Add Di Wu and Gordon Smyth as contributors, as much of the core functionality originates from their default method.
- Update `cluster_sets` function to return a `data.frame`, rather than a `data.table`.


# ostRich 0.1.0 (2024-04-22)

- Initial release.
