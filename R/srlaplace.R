

#' Sampling from Laplace Distribution
#' @rdname srlaplace
#' @order 1
#'
#' @description
#' Sampling from the Laplace distribution using Stors.
#'
#' @details
#' The function \code{srlaplace()} is used for sampling from a standard Laplace distribution (location = 0, scale = 1).
#'
#' Laplace distribution has the density:
#' \deqn{ f(x | \mu, b) = \frac{1}{2b} \exp\left(-\frac{|x - \mu|}{b}\right) }
#' Where \eqn{\mu} is the location parameter and \eqn{b} is the scale parameter.
#'
#' @param n Integer sample size.
#'
#' @return
#' \code{srlaplace()} returns a sample of size \code{n} from a standard Laplace distribution.
#' @examples
#' # Optimize the grid for the Laplace distribution
#' grid_optimizer("srlaplace")
#' # Generating Samples from a Standard Laplace Distribution
#' # This example illustrates how to generate 10 samples from a standard Laplace distribution.
#' samples <- srlaplace(10)
#' print(samples)
#'
#'
#' @export
srlaplace <- function(n = 1, mu = 0, b = 1, x = NULL) {
  .Call(C_srlaplace_scaled_check, n, c(mu, b), x)
}


#' Sampling from Custom Laplace Distribution
#' @rdname srlaplace
#' @order 2
#'
#' @description
#' Sampling from any Laplace distribution specifying the location and scale.
#'
#' @param mu Scalar location parameter.
#' @param b Scalar scale parameter.
#'
#' @return
#' \code{srlaplace_scaled()} returns a sample of size \code{n} from a Laplace distribution with location \eqn{\mu} and scale \eqn{b}.
#'
#' @examples
#' # This example demonstrates how to generate 10,
#' # Generating Samples from a Laplace Distribution with Specific Location and Scale
#' # samples from a Laplace distribution with location = 2 and scale = 3.
#' samples <- srlaplace_scaled(n = 10, mu = 2, b = 3)
#' print(samples)
#'
#' @export
srlaplace_custom <- function(n = 1, x = NULL) {
  .Call(C_srlaplace_custom_check, n, x)
}

#' @export
srlaplace_optimize <- function(
    mu = NULL,
    b = NULL,
    xl = NULL,
    xr = NULL,
    steps = 4091,
    grid_range = NULL,
    theta = NULL,
    target_sample_size = 1000,
    verbose = FALSE,
    symmetric = NULL
) {

  dist_name <- "srlaplace"

  dendata <- pbgrids[[dist_name]]

  f_params <- list(mu = mu, b = b)

  if (dendata$scalable) {

    isnull <- sapply(f_params, is.null)

    if (all(isnull)) {
      cnum <- dendata$c_num
      grid_type <- "scaled"
    } else {
      cnum <- dendata$c_num + 1
      grid_type <- "custom"
    }

    f_params <- ifelse(isnull, dendata$std_params, f_params)

  } else {

    cnum <- dendata$c_num + 1
    grid_type <- "custom"
  }

  modes <- dendata$set_modes(f_params$mu)

  f <- dendata$create_f(f_params$mu, f_params$b)

  check_grid_opt_criteria(symmetric, cnum, dendata)

  grid_optimizer(dendata, dist_name, xl, xr, f, modes, f_params, steps,
                 grid_range, theta, target_sample_size,
                 grid_type, symmetric, cnum, verbose)

}
