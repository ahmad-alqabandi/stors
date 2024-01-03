
#' @export
laplace <- function(){
  
  .Call(C_laplace, n)

}