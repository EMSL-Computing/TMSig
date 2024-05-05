
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MOST: MOlecular Signatures Toolkit

<!-- badges: start -->

![R package
version](https://img.shields.io/github/r-package/v/EMSL-Computing/most?label=R%20package)
<!-- badges: end -->

The `most` R package contains tools to prepare and analyze *a priori*
molecular signatures, such as gene sets. Examples of molecular
signatures: Reactome pathways, Gene Ontology gene sets, phosphorylation
sites grouped by their known kinases, and metabolites/lipids grouped by
chemical subclasses (e.g., acyl carnitines, fatty acids). **In general,
functions in this package work with any named list of character
vectors.**

Below is a overview of some core functions:

- `gmt_to_list`: create a named list of sets from a GMT file.

- `incidence`: compute a sparse incidence matrix with unique sets as
  rows and unique elements as columns. A value of 1 indicates that a
  particular element is a member of that set, while a value of 0
  indicates that it is not.

- `similarity`: compute the matrix of pairwise Jaccard or overlap
  similarity coefficients for all pairs of sets $A$ and $B$, where

  - $Jaccard(A, B) = \frac{|A \cap B|}{|A \cup B|}$
  - $Overlap(A, B) = \frac{|A \cap B|}{\min(|A|, |B|)}$

- `filter_sets`: restrict sets to only those elements in a
  pre-determined background, if provided, and only keep those that pass
  minimum and maximum size thresholds.

- `cluster_sets`: hierarchical clustering of highly similar sets. Used
  to reduce redundancy prior to analysis.

- `cameraPR.matrix`: a matrix method for `limma::cameraPR` for testing
  molecular signatures across one or more contrasts. Pre-Ranked
  Correlation Adjusted MEan RAnk molecular signature analysis
  (CAMERA-PR) accounts for inter-molecular correlation to control the
  type I error rate (Wu & Smyth, 2012).

- `enrichmap`: visualize molecular signature analysis results, such as
  those from `cameraPR.matrix`, as a bubble heatmap with signatures as
  rows and contrasts as columns.

## Installation

You can install the development version of `most` like so:

``` r
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("EMSL-Computing/most")
```

## Examples

``` r
library(most)
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
## Calculate matrices of pairwise Jaccard and overlap similarity coefficients
similarity(x) # Jaccard (default)
#> 7 x 7 sparse Matrix of class "dgCMatrix"
#>      Set1      Set2      Set3      Set4      Set5 Set6      Set7
#> Set1  1.0 0.8000000 0.8000000 0.6000000 0.2000000    . 0.5000000
#> Set2  0.8 1.0000000 1.0000000 0.7500000 0.2500000    . 0.3333333
#> Set3  0.8 1.0000000 1.0000000 0.7500000 0.2500000    . 0.3333333
#> Set4  0.6 0.7500000 0.7500000 1.0000000 0.3333333    . 0.1666667
#> Set5  0.2 0.2500000 0.2500000 0.3333333 1.0000000    . .        
#> Set6  .   .         .         .         .            1 .        
#> Set7  0.5 0.3333333 0.3333333 0.1666667 .            . 1.0000000

similarity(x, type = "overlap") # overlap
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
## Cluster sets based on their similarity

# Cluster aliased sets (the empty brackets are required to display the table)
cluster_sets(x, cutoff = 1)[]
#>       set cluster set_size
#>    <char>   <int>    <int>
#> 1:   Set2       1        4
#> 2:   Set3       1        4
#> 3:   Set1       2        5
#> 4:   Set4       3        3
#> 5:   Set5       4        1
#> 6:   Set6       5        3
#> 7:   Set7       6        4

# Cluster subsets
cluster_sets(x, cutoff = 1, type = "overlap")[]
#>       set cluster set_size
#>    <char>   <int>    <int>
#> 1:   Set1       1        5
#> 2:   Set2       1        4
#> 3:   Set3       1        4
#> 4:   Set4       1        3
#> 5:   Set5       1        1
#> 6:   Set6       2        3
#> 7:   Set7       3        4
```

## References

Wu, D., and Smyth, G. K. (2012). Camera: a competitive gene set test
accounting for inter-gene correlation. 40, e133. .
