#' @title Cluster similar sets
#'
#' @description Determine clusters of highly similar sets. Useful to reduce the
#'   redundancy of sets prior to statistical analysis.
#'
#' @param x a named list of sets. Only sets with at least 2 elements will be
#'   kept.
#' @inheritParams similarity
#' @param cutoff numeric 0-1; minimum value used to consider two sets as being
#'   similar.
#' @param h numeric 0-1; cut height used to define clusters. Passed to
#'   \code{\link[stats]{cutree}}. Default is 0.9, which seems to work reasonably
#'   well.
#'
#' @section Function Details:
#'
#'   Given a named list of sets, \code{cluster_sets} calculates all pairwise
#'   Jaccard or overlap similarity coefficients (see \code{\link{similarity}}
#'   for details). Any coefficients below \code{cutoff} are set to 0 and
#'   average-linkage hierarchical clustering is performed on the dissimilarity
#'   matrix (calculated as 1 - coefficients). Lastly,
#'   \code{\link[stats]{cutree}} is used with cut height \code{h} to define
#'   clusters, and the results are stored in a \code{data.frame}.
#'
#'   Results are arranged in ascending order by cluster, descending order by set
#'   size, and then alphanumerically by set name.
#'
#' @section Set Size:
#'
#'   Sets are pre-filtered to have at least 2 elements with
#'   \code{\link{filter_sets}}.
#'
#'   For two non-aliased sets to have \eqn{Jaccard \geq x}, where \eqn{x} is
#'   some threshold, sets must have minimum sizes \eqn{n1 = \lceil \frac{x}{1 -
#'   x} \rceil} and \eqn{n2 = n1 + 1} with the size of the overlap equal to
#'   \eqn{n1}. For example, if sets have \eqn{J \geq 0.85}, their minimum sizes
#'   are 6 and 7 with an overlap size of 6. That is, each set with fewer than 6
#'   elements will always appear as a singleton cluster.
#'
#'   For two sets to have \eqn{Overlap \geq x}, where one set is not a subset of
#'   the other and \eqn{x} is some threshold, the sets must both be of minimum
#'   size \eqn{n = 1 + \lceil \frac{x}{1 - x} \rceil} with the size of the
#'   overlap equal to \eqn{n - 1}. For example, if sets have \eqn{Overlap \geq
#'   0.85}, their minimum sizes are both 7 with an overlap size of 6. That is,
#'   each set with fewer than 7 elements will always appear as a singleton
#'   cluster.
#'
#' @section Optimization:
#'
#'   Clustering does not need to be performed on those sets that are not
#'   sufficiently similar (value of \code{similarity} below \code{cutoff}) to
#'   any other set, as they will always be placed in their own cluster. By
#'   excluding these sets during the hierarchical clustering step, the speed of
#'   \code{cluster_sets} will increase as the value of \code{cutoff} approaches
#'   1 (as the size of the dissimilarity matrix decreases).
#'
#' @source
#'
#' This function is based on the procedure described in the Molecular Signatures
#' Database (MSigDB) v7.0 Release Notes (Liberzon 2011, 2015):
#' \url{https://docs.gsea-msigdb.org/#MSigDB/Release_Notes/MSigDB_7.0/}.
#' Specifically, sections "C2:CP:Reactome — Major Overhaul" and "C5 (Gene
#' Ontology Collection) — Major Overhaul". Though defining clusters from a
#' dissimilarity matrix is a commonly used approach, the default values of
#' \code{cutoff} and cut height \code{h} are exactly what is set for MSigDB
#' (private correspondence).
#'
#' @return A \code{data.frame} with columns \code{"set"}, \code{"cluster"}, and
#'   \code{"set_size"}.
#'
#' @references
#'
#' Liberzon, A., Subramanian, A., Pinchback, R., Thorvaldsdóttir, H., Tamayo,
#' P., & Mesirov, J. P. (2011). Molecular signatures database (MSigDB) 3.0.
#' \emph{Bioinformatics, 27}(12), 1739–1740.
#' \url{https://doi.org/10.1093/bioinformatics/btr260}
#'
#' Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J. P., &
#' Tamayo, P. (2015). The Molecular Signatures Database (MSigDB) hallmark gene
#' set collection. \emph{Cell systems, 1}(6), 417–425.
#' \url{https://doi.org/10.1016/j.cels.2015.12.004}
#'
#' @seealso \code{\link{filter_sets}}, \code{\link{similarity}}
#'
#' @examples
#' x <- list("A" = letters[1:5],
#'           "B" = letters[1:4], # subset of A
#'           "C" = letters[1:4], # aliased with B
#'           "D" = letters[1:3], # subset of A, B, C
#'           "E" = c("a", "a", NA), # duplicates and NA
#'           "F" = c("x", "y", "z"), # distinct elements
#'           "G" = letters[3:6]) # overlaps with A-E
#'
#' cluster_sets(x)
#'
#' # Relax similarity cutoff
#' (df <- cluster_sets(x, cutoff = 0.5))
#'
#' # Keep the first (largest) set from each cluster
#' keep <- df[!duplicated(df$cluster), "set"]
#' x[keep] # A, F
#'
#' # Keep the smallest set from each cluster
#' df <- df[order(df$set_size), ] # preserves order of set names
#' keep <- df[!duplicated(df$cluster), "set"]
#' x[keep] # D, F
#'
#' # Cluster aliased sets
#' cluster_sets(x, method = "jaccard", cutoff = 1)
#'
#' # Cluster subsets
#' cluster_sets(x, method = "overlap", cutoff = 1)
#'
#' @importFrom stats as.dist hclust cutree
#'
#' @export cluster_sets

cluster_sets <- function(x,
                         method = c("jaccard", "overlap"),
                         cutoff = 0.85,
                         h = 0.9)
{
  if (cutoff < 0 | cutoff > 1)
    stop("`cutoff` must be between 0 and 1.")

  # Pre-filter sets (also validates x)
  x <- filter_sets(x, min_size = 2L)

  if (length(x) < 2L)
    stop("`x` must contain 2 or more sets.")

  # Similarity matrix
  s <- similarity(x, method = method)
  diag(s) <- 0
  s[s < cutoff] <- 0

  # Cluster sets that are similar to at least one other set. The remaining sets
  # will be appended to the results at the end.
  keep <- apply(s, 1, function(.s) any(.s > 0))

  if (sum(keep) == 0L) {
    message("No pair of sets passes the similarity cutoff.")

    df <- data.frame(set = names(x),
                     cluster = seq_along(x),
                     set_size = lengths(x),
                     row.names = NULL,
                     stringsAsFactors = FALSE)

    o <- with(df, order(cluster, set_size, set, method = "radix",
                        decreasing = c(FALSE, TRUE, FALSE)))
    df <- df[o, ]
    rownames(df) <- NULL # order() modifies rownames

    return(df)
  }

  s <- s[keep, keep] # at least a 2x2 matrix

  # Convert sparse similarity matrix to dense dissimilarity matrix
  d <- (1 - s)
  d <- as.dist(d)

  # Hierarchical clustering
  hc <- hclust(d, method = "average")
  clusters <- cutree(hc, h = h)

  df <- data.frame(set = names(clusters),
                   cluster = clusters,
                   row.names = NULL,
                   stringsAsFactors = FALSE)

  # Append sets that were not similar to any others and place each in its own
  # cluster
  other_sets <- data.frame(set = names(x)[!keep],
                           row.names = NULL,
                           stringsAsFactors = FALSE)
  other_sets$cluster <- max(df$cluster) + seq_len(nrow(other_sets))

  df <- rbind(df, other_sets)

  df$set_size <- lengths(x)[df$set]

  o <- with(df, order(cluster, set_size, set, method = "radix",
                      decreasing = c(FALSE, TRUE, FALSE)))
  df <- df[o, ]
  rownames(df) <- NULL # order() modifies rownames

  return(df)
}

