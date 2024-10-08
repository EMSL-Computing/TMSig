# TMSig DEV (2024-09-07)

**SIGNIFICANT USER-VISIBLE CHANGES**

- Convert functions `sparseIncidence`, `filterSets`, and `invertSets` to S4 generics and add methods for objects of class `GSEABase::GeneSet` and `GSEABase::GeneSetCollection`.
- Rename functions to adhere to Bioconductor style: `cluster_sets` to `clusterSets`, `decompose_sets` to `decomposeSets`, `filter_sets` to `filterSets`, `gmt_to_list` to `readGMT`, `incidence` to `sparseIncidence`, `incidence_to_list` to `incidenceToList`, `invert_sets` to `invertSets`, `range_extend` to `extendRangeNum`, `camera_color_fun` to `cameraColorFun`, `gsea_color_fun` to `gseaColorFun`.
- `enrichmap()`: heatmap rectangle fill and border color can now be changed via the `heatmap_args` parameter. 
- `enrichmap()`: color legend title can now be changed via `heatmap_args = list(name = "new title")` instead of `heatmap_args = list(heatmap_legend_param = list(title = "new title"))`. 
- `enrichmap()`: added examples to documentation. Default text size is now controlled by the `cell_size` argument.
- `decomposeSets()`: added arguments `AND`, `MINUS`, and `verbose`. Reduced computation time.
- "GO_analysis" vignette replaced by "TMSig" vignette.

**MINOR CHANGES**

- Simplified internals of several functions. In most cases, this also led to a slight reduction in computation time and memory usage.
- Minor documentation fixes.

**BUG FIXES**

- Added Bioconductor packages to Remotes field of DESCRIPTION.
- `cameraPR.matrix()`: fixed integer overflow that could occur when testing very large sets with `use.ranks=TRUE`. 
- `cameraPR.matrix()`: fixed negative variance caused by improper adjustment for ties in ranks when `use.ranks=TRUE`.


# TMSig 1.1.1 (2024-07-02)

**BUG FIXES**

- `cameraPR.matrix()`: arguments `before` and `after` of `data.table::setcolorder()` were added in data.table v1.15.0 (2024-01-30). Rather than requiring a recent minimum version of data.table, the code has been modified to no longer depend on these arguments.


# TMSig 1.1.0 (2024-06-26)

**SIGNIFICANT USER-VISIBLE CHANGES**

- `cameraPR.matrix()`: added `alternative` parameter to allow users to perform one-sided tests. A warning will be issued if attempting to use a one-sided test with the parametric version of CAMERA-PR, as this is generally not recommended.
- `cameraPR.matrix()`: added `min.size` parameter to allow users to specify a minimum set size for testing.
- `similarity()` and `cluster_sets()`: added support for Ōtsuka set similarity.
- Documentation for most functions has been updated.

**BUG FIXES**

- `cameraPR.matrix()`: fixed variance calculation when `use.ranks=TRUE` and `inter.gene.cor=0`.
- `enrichmap()`: fixed checks for duplicate contrast-set pairs. 
- `enrichmap()`: if column or row labels are specified by the user through the `heatmap_args` parameter, the maximum text width or height is now determined from these new labels to avoid overlapping plot elements.
- `enrichmap`: default of `statistic_column` parameter changed to "TwoSampleT" to match output of `cameraPR.matrix`.


# TMSig 1.0.0 (2024-05-06)

**SIGNIFICANT USER-VISIBLE CHANGES**

- Rename package from ostRich to TMSig: Tools for Molecular Signatures.
- Rename "prepare_gene_sets" vignette to "GO_analysis". Add example code for `cameraPR.matrix`.

**NEW FEATURES**

- Add fast `cameraPR.matrix` method based on `limma::cameraPR.default`. Add Di Wu and Gordon Smyth as contributors, as much of the core functionality originates from their default method.
- Update `cluster_sets` function to return a `data.frame`, rather than a `data.table`.


# ostRich 0.1.0 (2024-04-22)

- Initial release.
