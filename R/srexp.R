#' @export
srexp <- function(n) {
  .Call(C_srexp, n)
}
