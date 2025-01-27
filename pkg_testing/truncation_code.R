

# srnorm_truncate <- function(xl = -Inf, xr = Inf, mean = 0, sd = 1){
#
#   dist_name <- 'srnorm'
#
#   dendata <- built_in_proposals[[dist_name]]
#
#   cnum_scalable <- dendata$c_num
#   cnum_custom <- dendata$c_num + 1
#   choosen_grid_num <- NULL
#
#   scalable_info <-  cached_proposal_info(cnum_scalable)
#   custom_info <-  cached_proposal_info(cnum_custom)
#
#   res <- truncate_error_checking(xl, xr, dendata)
#   xl <- res$xl; xr <- res$xr
#
#   if( !is.null(custom_info) && all(custom_info[-1] == c(mean, sd)) ){
#     choosen_grid_num <- cnum_custom
#     Upper_cumsum <- .Call(C_srnorm_trunc_nav, xl, xr, choosen_grid_num)
#
#   }else if( !is.null(scalable_info)){
#     choosen_grid_num <- cnum_scalable
#     Upper_cumsum <- .Call(C_srnorm_trunc_nav, xl, xr, choosen_grid_num)
#
#   } else{
#     .Call(C_grid_error,0,0)
#     return(NULL)
#
#   }
#
#
#   stopifnot(
#     "xl is has a CDF close to 1" = (Upper_cumsum[1] != 1),
#     "xr is has a CDF close to 0" = (Upper_cumsum[2] != 0)
#   )
#
#   function_string <- paste0("function(n) { .Call(C_srnorm_trunc, n, ",
#                             paste0(xl), ", ", paste0(xr), ", ", paste0(Upper_cumsum[1]),
#                             ", ", paste0(Upper_cumsum[2]), ", ", paste0(as.integer(Upper_cumsum[3])), ", ",
#                             paste0(as.integer(Upper_cumsum[4])),", ",
#                             paste0(as.integer(choosen_grid_num)), ")}")
#
#   function_expression <- parse(text = function_string)
#   sampling_function <- eval(function_expression)
#
#   return(sampling_function)
# }

