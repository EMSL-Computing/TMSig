#' @title Similarity-Based Clustering of Sets
#'
#' @description Determine clusters of highly similar sets. Used to reduce the
#'   redundancy of sets prior to statistical analysis.
#'
#' @inheritParams similarity
#' @param cutoff numeric 0-1; minimum similarity coefficient required to
#'   classify two sets as being similar. Default is 0.85.
#' @param method character; the clustering method passed to
#'   \code{\link[stats]{hclust}}. Default is "complete", so sets will only be
#'   included in a cluster if their similarity to all other sets in that cluster
#'   is \eqn{\geq} \code{cutoff}.
#' @param h numeric 0-1; cut height used to define clusters. Passed to
#'   \code{\link[stats]{cutree}}. Default is 0.9.
#'
#' @returns A \code{data.frame} with 3 columns:
#'
#'   \item{\code{set}}{character; the name of the set.}
#'   \item{\code{cluster}}{integer; the cluster identifier.}
#'   \item{\code{set_size}}{integer; the size of the set (number of elements).}
#'
#'   Results are arranged in ascending order by cluster, descending order by set
#'   size, and then alphanumerically by set name.
#'
#' @section Function Details:
#'
#'   Given a named list of sets, \code{cluster_sets} calculates all pairwise
#'   Jaccard or overlap similarity coefficients (see \code{\link{similarity}}
#'   for details). Any coefficients below \code{cutoff} are set to 0 and
#'   complete-linkage hierarchical clustering is performed on the dissimilarity
#'   matrix (calculated as 1 - coefficients). Lastly,
#'   \code{\link[stats]{cutree}} is used with cut height \code{h} to define
#'   clusters, and the results are stored in a \code{data.table}.
#'
#' @section Set Size:
#'
#'   For two non-aliased sets to have \eqn{Jaccard \geq x}, where \eqn{x} is
#'   some threshold, sets must have minimum sizes \eqn{n1 = \lceil \frac{x}{1 -
#'   x} \rceil} and \eqn{n2 = n1 + 1} with the size of the overlap equal to
#'   \eqn{n1}. For example, if sets have \eqn{J \geq 0.85}, their minimum sizes
#'   are 6 and 7 with an overlap size of 6. That is, each set with fewer than 6
#'   elements will always appear as a singleton cluster, if they are not
#'   aliased.
#'
#'   For two non-aliased sets to have \eqn{Overlap \geq x}, where one set is not
#'   a subset of the other and \eqn{x} is some threshold, the sets must both be
#'   of minimum size \eqn{n = 1 + \lceil \frac{x}{1 - x} \rceil} with the size
#'   of the overlap equal to \eqn{n - 1}. For example, if sets have \eqn{Overlap
#'   \geq 0.85}, their minimum sizes are both 7 with an overlap size of 6. That
#'   is, each set with fewer than 7 elements will always appear as a singleton
#'   cluster, if they are not subsets or aliased.
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
#' Ontology Collection) — Major Overhaul". Though hierarchical clustering is
#' widely used, the defaults are exactly what is set for MSigDB (private
#' correspondence).
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
#' @importFrom data.table data.table := setorderv setDF
#' @importFrom stats as.dist hclust cutree
#'
#' @export cluster_sets
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
#' cluster_sets(x) # default clustering
#'
#' # Relax similarity cutoff
#' (df <- cluster_sets(x, cutoff = 0.5))
#'
#' # Keep the first (largest) set from each cluster
#' subset(df, subset = !duplicated(cluster))[["set"]] # A, G, E, F
#'
#' # Keep the smallest set from each cluster
#' df <- df[order(df$set_size), ]
#' subset(df, subset = !duplicated(cluster))[["set"]] # E, D, F, G
#'
#' # Cluster aliased sets
#' cluster_sets(x, type = "jaccard", cutoff = 1)
#'
#' # Cluster subsets
#' cluster_sets(x, type = "overlap", cutoff = 1)

cluster_sets <- function(x,
                         type = c("jaccard", "overlap"),
                         cutoff = 0.85,
                         method = "complete",
                         h = 0.9)
{
  if (cutoff < 0 | cutoff > 1)
    stop("`cutoff` must be between 0 and 1.")

  # x is used later to get set sizes
  dt <- .prepare_sets(x)
  x <- split(x = dt[["elements"]], f = dt[["sets"]])

  # Similarity matrix (not optimal because similarity uses .prepare_sets)
  s <- similarity(x, type = type)
  diag(s) <- 0
  s[s < cutoff] <- 0

  # Cluster sets that are sufficiently similar to at least one other set. The
  # remaining sets will be appended to the results at the end.
  keep <- apply(s > 0, 1, any)

  if (sum(keep) == 0L) {
    message("No pair of sets passes the similarity cutoff.")

    dt <- data.table(set = names(x),
                     cluster = seq_along(x),
                     set_size = lengths(x),
                     stringsAsFactors = FALSE)

    setorderv(dt, cols = c("cluster", "set_size", "set"),
              order = c(1, -1, 1))

    setDF(dt) # convert to data.frame

    return(dt)
  }

  s <- s[keep, keep] # at least a 2x2 matrix

  # Convert sparse similarity matrix to dense dissimilarity matrix - may produce
  # a warning
  d <- as.dist(1 - s)

  # Hierarchical clustering
  hc <- hclust(d, method = method)
  clusters <- cutree(hc, h = h)

  dt <- data.table(set = names(clusters),
                   cluster = clusters,
                   stringsAsFactors = FALSE)

  # Append sets that were not similar to any others and place each in its own
  # cluster
  other_sets <- data.table(set = names(x)[!keep],
                           stringsAsFactors = FALSE)
  other_sets[, cluster := seq_along(set) + max(dt$cluster)]

  dt <- rbind(dt, other_sets)
  dt[, set_size := lengths(x)[set]]

  setorderv(dt, cols = c("cluster", "set_size", "set"),
            order = c(1, -1, 1))

  setDF(dt) # convert to data.frame

  return(dt)
}
