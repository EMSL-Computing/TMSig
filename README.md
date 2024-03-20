
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ostRich

<!-- badges: start -->

![R package
version](https://img.shields.io/github/r-package/v/EMSL-Computing/ostRich?label=R%20package)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

ostRich contains tools to operate on sets, and it was originally
intended for the preparation of biomolecular/-omics sets for enrichment
analysis/gene set testing/biomolecular set testing. Examples of
biomolecular sets: Reactome pathways, Gene Ontology gene sets,
phosphorylation sites grouped by their known kinases, and
metabolites/lipids grouped by chemical subclasses (e.g., acyl
carnitines, fatty acids). **In general, the functions in this package
work with any named list of character vectors.**

Below is a growing list of what can be accomplished with this package:

- `gmt_to_list`: Create a named list of gene sets from a GMT file.

- `incidence`: Compute a sparse incidence matrix with unique sets as
  rows and unique elements as columns. A value of 1 indicates that a
  particular element is a member of that set, while a value of 0
  indicates that it is not. This serves as the cornerstone of several
  other functions. The cross-product of the incidence matrix with its
  transpose will yield the sizes of the pairwise intersections of each
  set on the off-diagonals and the sizes of the sets on the diagonal.

- `similarity`: Compute all pairwise Jaccard or overlap similarity
  coefficients. The Jaccard coefficient, $J$, is defined as the
  cardinality of the intersection of two sets divided by the cardinality
  of their union. The overlap coefficient is defined as the cardinality
  of the intersection divided by the cardinality of the smallest set.

- `filter_sets`: optionally restrict sets to only those elements in a
  pre-determined `background` and only keep those sets passing minimum
  and maximum size thresholds. In practice, biomolecular sets are
  restricted to the biomolecules measured in a particular experiment
  (the background).

- `cluster_sets`: hierarchical clustering of highly similar sets. In
  practice, a single biomolecular set from each cluster would be
  selected to serve as that cluster’s representative for enrichment
  analysis. Redundant sets would be discarded or reported separately
  (preferred).

- `decompose_sets`: decompose all pairs of sets into 3 disjoint
  components: the elements unique to set 1, the elements unique to set
  2, and the elements shared by both sets. Not currently used in
  practice.

## Installation

You can install the development version of ostRich like so:

``` r
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("EMSL-Computing/ostRich")
```

## Examples

``` r
library(ostRich)
#> Loading required package: Matrix

# Named list of sets
x <- list("Set1" = letters[1:5],
          "Set2" = letters[1:4], # subset of Set1
          "Set3" = letters[1:4], # aliased with Set2
          "Set4" = letters[1:3], # subset of Set1-Set3
          "Set5" = c("a", "a", NA), # duplicates and NA
          "Set6" = c("x", "y", "z"), # distinct elements
          "Set7" = letters[3:6]) # overlaps with Set1-Set5
x
#> $Set1
#> [1] "a" "b" "c" "d" "e"
#> 
#> $Set2
#> [1] "a" "b" "c" "d"
#> 
#> $Set3
#> [1] "a" "b" "c" "d"
#> 
#> $Set4
#> [1] "a" "b" "c"
#> 
#> $Set5
#> [1] "a" "a" NA 
#> 
#> $Set6
#> [1] "x" "y" "z"
#> 
#> $Set7
#> [1] "c" "d" "e" "f"
```

``` r
(imat <- incidence(x)) # incidence matrix
#> 7 x 9 sparse Matrix of class "dgCMatrix"
#>      a b c d e x y z f
#> Set1 1 1 1 1 1 . . . .
#> Set2 1 1 1 1 . . . . .
#> Set3 1 1 1 1 . . . . .
#> Set4 1 1 1 . . . . . .
#> Set5 1 . . . . . . . .
#> Set6 . . . . . 1 1 1 .
#> Set7 . . 1 1 1 . . . 1

tcrossprod(imat) # pairwise intersection and set sizes
#> 7 x 7 sparse Matrix of class "dsCMatrix"
#>      Set1 Set2 Set3 Set4 Set5 Set6 Set7
#> Set1    5    4    4    3    1    .    3
#> Set2    4    4    4    3    1    .    2
#> Set3    4    4    4    3    1    .    2
#> Set4    3    3    3    3    1    .    1
#> Set5    1    1    1    1    1    .    .
#> Set6    .    .    .    .    .    3    .
#> Set7    3    2    2    1    .    .    4

crossprod(imat) # occurrence of each element and pair of elements
#> 9 x 9 sparse Matrix of class "dsCMatrix"
#>   a b c d e x y z f
#> a 5 4 4 3 1 . . . .
#> b 4 4 4 3 1 . . . .
#> c 4 4 5 4 2 . . . 1
#> d 3 3 4 4 2 . . . 1
#> e 1 1 2 2 2 . . . 1
#> x . . . . . 1 1 1 .
#> y . . . . . 1 1 1 .
#> z . . . . . 1 1 1 .
#> f . . 1 1 1 . . . 1
```

``` r
similarity(x) # pairwise Jaccard similarity
#> 7 x 7 sparse Matrix of class "dgCMatrix"
#>      Set1      Set2      Set3      Set4      Set5 Set6      Set7
#> Set1  1.0 0.8000000 0.8000000 0.6000000 0.2000000    . 0.5000000
#> Set2  0.8 1.0000000 1.0000000 0.7500000 0.2500000    . 0.3333333
#> Set3  0.8 1.0000000 1.0000000 0.7500000 0.2500000    . 0.3333333
#> Set4  0.6 0.7500000 0.7500000 1.0000000 0.3333333    . 0.1666667
#> Set5  0.2 0.2500000 0.2500000 0.3333333 1.0000000    . .        
#> Set6  .   .         .         .         .            1 .        
#> Set7  0.5 0.3333333 0.3333333 0.1666667 .            . 1.0000000

similarity(x, method = "overlap") # pairwise overlap similarity
#> 7 x 7 sparse Matrix of class "dgCMatrix"
#>      Set1 Set2 Set3      Set4 Set5 Set6      Set7
#> Set1 1.00  1.0  1.0 1.0000000    1    . 0.7500000
#> Set2 1.00  1.0  1.0 1.0000000    1    . 0.5000000
#> Set3 1.00  1.0  1.0 1.0000000    1    . 0.5000000
#> Set4 1.00  1.0  1.0 1.0000000    1    . 0.3333333
#> Set5 1.00  1.0  1.0 1.0000000    1    . .        
#> Set6 .     .    .   .            .    1 .        
#> Set7 0.75  0.5  0.5 0.3333333    .    . 1.0000000
```

``` r
cluster_sets(x, cutoff = 1) # cluster aliased sets
#>    set cluster set_size
#> 1 Set2       1        4
#> 2 Set3       1        4
#> 3 Set1       2        5
#> 4 Set4       3        3
#> 5 Set6       4        3
#> 6 Set7       5        4

cluster_sets(x, cutoff = 1, method = "overlap") # cluster subsets
#>    set cluster set_size
#> 1 Set1       1        5
#> 2 Set2       1        4
#> 3 Set3       1        4
#> 4 Set4       1        3
#> 5 Set6       2        3
#> 6 Set7       3        4
```
