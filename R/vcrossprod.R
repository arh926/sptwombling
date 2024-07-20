#' Fucntion for vector cross-product
#'
#' Computes a cross product for two vectors
#'
#' @param u vector 1
#' @param v vector 2
#' @keywords
#' @export
#' @examples
#################################
# Cross-product for two vectors #
#################################
vcrossprod <- function(u = NULL,
                       v = NULL){
  result = c(u[2] * v[3] - u[3] * v[2],
             -(u[1] * v[3] - u[3] * v[1]),
             u[1] * v[2] - u[2] * v[1])
  names(result) = c("x", "y", "t")
  return(result)
}
