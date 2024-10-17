# TMSig 0.99.2 (2024-10-17)

**MINOR CHANGES**
- Reduced size of TMSig vignette.


# TMSig 0.99.1 (2024-09-30)

**SIGNIFICANT USER-VISIBLE CHANGES**

- Renamed functions to adhere to Bioconductor style: `cluster_sets` to `clusterSets`, `decompose_sets` to `decomposeSets`, `filter_sets` to `filterSets`, `gmt_to_list` to `readGMT`, `incidence` to `sparseIncidence`, `incidence_to_list` to `incidenceToList`, `invert_sets` to `invertSets`, `range_extend` to `extendRangeNum`, `camera_color_fun` to `cameraColorFun`, `gsea_color_fun` to `gseaColorFun`.
- Converted functions `sparseIncidence`, `filterSets`, and `invertSets` to S4 generics and add methods for objects of class `GSEABase::GeneSet` and `GSEABase::GeneSetCollection`. Added paragraph explaining how `sparseIncidence` differs from `GSEABase::incidence`.
- `cameraPR.matrix()`: added "ZScore" column to results (placed before "PValue" column). The default value of the `statistic_column` parameter of `enrichmap` was also updated to `"ZScore"`.
- `enrichmap()`: heatmap rectangle fill and border color can now be changed via the `heatmap_args` parameter.
- `enrichmap()`: color legend title can now be changed via `heatmap_args = list(name = "new title")` instead of `heatmap_args = list(heatmap_legend_param = list(title = "new title"))`. 
- `enrichmap()`: added examples to documentation. Default text size is now controlled by the `cell_size` argument.
- `enrichmap()`: updated the heatmap layer function; all bubbles, regardless of significance, will have a radius of at least `0.2 * cell_size`.

**MINOR CHANGES**

- Added inst/script/c5.go.v2023.2.Hs.symbols.txt file to describe the process used to generate c5.go.v2023.2.Hs.symbols.gmt.gz.

**BUG FIXES**

- Fixed manual build errors caused by special characters used in LaTeX equations.


# TMSig 0.99.0 (2024-08-24)

- Initial version.
