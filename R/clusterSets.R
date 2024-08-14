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
#'   Given a named list of sets, \code{clusterSets} calculates all pairwise
#'   Jaccard, overlap, or Ōtsuka similarity coefficients (see
#'   \code{\link{similarity}} for details). Any coefficients below \code{cutoff}
#'   are set to 0 and complete-linkage hierarchical clustering is performed on
#'   the dissimilarity matrix (calculated as 1 - coefficients). Lastly,
#'   \code{\link[stats]{cutree}} is used with cut height \code{h} to define
#'   clusters, and the results are stored in a \code{data.frame}.
#'
#' @section Optimization:
#'
#'   Clustering does not need to be performed on those sets that are not
#'   sufficiently similar (value of \code{similarity} below \code{cutoff}) to
#'   any other set, as they will always be placed in their own cluster. By
#'   excluding these sets during the hierarchical clustering step, the speed of
#'   \code{clusterSets} will increase as the value of \code{cutoff} approaches
#'   1 (as the size of the dissimilarity matrix decreases).
#'
#' @section Minimum Set Size:
#'
#'   Sets that are not sufficiently large will always appear as singleton
#'   clusters, unless they are aliased or subsets (overlap similarity only). For
#'   two sets \eqn{A} and \eqn{B} to be sufficiently similar, defined as having
#'   a similarity coefficient at least equal to some cutoff (e.g.,
#'   \eqn{Jaccard~\geq~x}), they must have minimum sizes \eqn{|A|}, \eqn{|B|},
#'   and intersection size \eqn{|A \cap B|}:
#'
#'   \itemize{
#'     \item{\strong{Jaccard:} \eqn{|A| = \lceil \frac{x}{1 - x} \rceil,
#'     \quad |B| = 1 + |A|, \quad |A \cap B| = |A|}}
#'     \item{\strong{Overlap:} \eqn{|A| = |B| = 1 + \lceil \frac{x}{1 - x}
#'     \rceil, \quad |A \cap B| = |A| - 1}}
#'     \item{\strong{Ōtsuka:} \eqn{|A| = \lceil \frac{x^2}{1 - x^2} \rceil,
#'     \quad |B| = 1 + |A|, \quad |A \cap B| = |A|}}
#'   }
#'
#'   where \eqn{\lceil y \rceil} is the ceiling function applied to some real
#'   number \eqn{y}.
#'
#'   For example, if the cutoff is \eqn{x = 0.85}, then the minimum set and
#'   intersection sizes are
#'
#'   \itemize{
#'     \item{\strong{Jaccard:} \eqn{|A| = 6, \quad |B| = 7, \quad |A \cap B| =
#'     6}}
#'     \item{\strong{Overlap:} \eqn{|A| = |B| = 7, \quad |A \cap B| = 6}}
#'     \item{\strong{Ōtsuka:} \eqn{|A| = 3, \quad |B| = 4, \quad |A \cap B| =
#'     3}}
#'   }
#'
#'   That is, sets with fewer elements or smaller intersections will always
#'   appear as singleton clusters unless they are aliased or, in the case of the
#'   overlap similarity, subsets.
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
#' doi:\href{https://doi.org/10.1093/bioinformatics/btr260
#' }{10.1093/bioinformatics/btr260}
#'
#' Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J. P., &
#' Tamayo, P. (2015). The Molecular Signatures Database (MSigDB) hallmark gene
#' set collection. \emph{Cell systems, 1}(6), 417–425.
#' doi:\href{https://doi.org/10.1016/j.cels.2015.12.004
#' }{10.1016/j.cels.2015.12.004}
#'
#' @seealso \code{\link{filterSets}}, \code{\link{similarity}}
#'
#' @importFrom data.table data.table := setorderv setDF
#' @importFrom Matrix diag
#' @importFrom stats as.dist hclust cutree
#'
#' @export clusterSets
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
#' # Default clustering based on Jaccard similarity
#' clusterSets(x)
#' clusterSets(x, type = "overlap") # overlap similarity
#' clusterSets(x, type = "otsuka") # Ōtsuka similarity
#'
#' # Relax Jaccard similarity cutoff
#' (df <- clusterSets(x, cutoff = 0.5))
#'
#' # Keep the first (largest) set from each cluster
#' with(df, set[!duplicated(cluster)]) # A, G, E, F
#'
#' # Keep the smallest set from each cluster
#' df <- df[order(df$set_size), ]
#' with(df, set[!duplicated(cluster)]) # E, D, F, G
#'
#' # Cluster aliased sets (type = "otsuka" would produce
#' # identical results)
#' clusterSets(x, type = "jaccard", cutoff = 1)
#'
#' # Cluster subsets and aliased sets
#' clusterSets(x, type = "overlap", cutoff = 1)

clusterSets <- function(x,
                        type = c("jaccard", "overlap", "otsuka"),
                        cutoff = 0.85,
                        method = "complete",
                        h = 0.9) {
    if (cutoff < 0 | cutoff > 1)
        stop("`cutoff` must be between 0 and 1.")

    # Prepare x to determine set sizes
    dt <- .prepare_sets(x)
    x <- split(x = dt[["elements"]], f = dt[["sets"]])
    set_sizes <- lengths(x)

    ## Pairwise set similarity matrix
    s <- similarity(x, type = type)
    diag(s) <- 0
    s[s < cutoff] <- 0

    # Cluster sets that are sufficiently similar to at least one other set
    keep <- apply(s != 0, 1, any)

    if (sum(keep) == 0L) {
        message("No pair of sets passes the similarity cutoff.")

        dt <- data.table(set = names(set_sizes),
                         cluster = seq_along(set_sizes),
                         stringsAsFactors = FALSE)
    } else {
        s <- s[keep, keep] # at least a 2x2 matrix

        # Convert sparse similarity matrix to dense dissimilarity matrix
        d <- as.dist(1 - s) # may produce a warning

        # Hierarchical clustering
        hc <- hclust(d, method = method)
        clusters <- cutree(hc, h = h)

        dt <- data.table(set = names(clusters),
                         cluster = clusters,
                         stringsAsFactors = FALSE)

        # Sets not similar to any others are placed in their own clusters
        other_sets <- data.table(set = names(set_sizes)[!keep],
                                 stringsAsFactors = FALSE)
        other_sets[, cluster := seq_along(set) + max(dt[["cluster"]])]

        dt <- rbind(dt, other_sets)
    }

    dt[, set_size := set_sizes[set]]

    setorderv(dt, cols = c("cluster", "set_size", "set"), order = c(1, -1, 1))

    setDF(dt)

    return(dt)
}
