# TMSig 1.1.0

**Improvements**

- `cameraPR.matrix`: added `alternative` parameter to allow users to perform one-sided tests. A warning will be issued if attempting to use a one-sided test with the parametric version of CAMERA-PR, as this is generally not recommended.
- `cameraPR.matrix`: added `min.size` parameter to allow users to specify a minimum set size for testing.
- `similarity` and `cluster_sets`: added support for Ōtsuka set similarity.
- Documentation for most functions has been updated.

**Bugfixes**

- `cameraPR.matrix`: fixed variance calculation when `use.ranks=TRUE` and `inter.gene.cor=0`.
- `enrichmap`: fixed checks for duplicate contrast-set pairs. 
- `enrichmap`: if column or row labels are specified by the user through the `heatmap_args` parameter, the maximum text width or height is now determined from these new labels to avoid overlapping plot elements.
- `enrichmap`: default of `statistic_column` parameter changed to "TwoSampleT" to match output of `cameraPR.matrix`.


# TMSig 1.0.0

- Rename package from ostRich to TMSig: Tools for Molecular Signatures.
- Add fast `cameraPR.matrix` method based on `limma::cameraPR.default`. Add Di Wu and Gordon Smyth as contributors.
- Update `cluster_sets` function to return a `data.frame`, rather than a `data.table`.
- Rename `prepare_gene_sets` vignette to `GO_analysis`. Add example code for `cameraPR.matrix`.


# ostRich 0.1.0

- Initial release.
