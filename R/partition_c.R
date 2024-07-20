#' Partitions a planar curve in space
#'
#' @param curve.list a list of planar curves
#' @param length.p length of partition
#' @keywords
#' @export
#' @examples
######################
# Partitions a curve #
######################
partition_curve <- function(curve.list = NULL, length.p = 10){
  nt = length(curve.list)

  curves.part = do.call(rbind, sapply(1:nt, function(x){
    id.part = round(seq(1, nrow(curve.list[[x]]), length.out = length.p))
    cbind(curve.list[[x]][id.part,], t = x)
  }, simplify = FALSE))

  colnames(curves.part) = c("long", "lat", "t")
  curves.part
}
