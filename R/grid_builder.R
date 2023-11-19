# pbgrids : 1 Normal

pbgrids <- list(
  snorm = c(
    Cnum = 0,
    f = function(x) {
      0.3989423 * exp(-0.5 * x * x)
    },
    h = function(x) {
      log(0.3989423) - (x * x) * (1 / 2)
    },
    h_prime = function(x) {
      -x
    },
    modes = 0,
    lb = -Inf,
    rb = Inf
  )
)


#' Title
#'
#' @param density_name 
#'
#' @return
#' @export
#'
#' @examples
grid_obtimizer <- function(density_name = c("snorm")) {
  density_name <- match.arg(density_name)

  dendata <- pbgrids[[density_name]]

  A <- seq(from = 0.005, to = 0.0002, length.out = 100)

  for (a in A) {
    opt_grid <- grid_builder(a = a, th = 0.6, mode = dendata$modes, f = dendata$f, h = dendata$h, h_prime = dendata$h_prime)

    if ((object.size(opt_grid) / 1024) > 128) break
  }
  
  opt_grid <- grid_builder(a = 0.00035, th = 0.6, mode = dendata$modes, f = dendata$f, h = dendata$h, h_prime = dendata$h_prime)
  
  
  cash_grid( dendata$Cnum ,opt_grid)
  
}

cash_grid = function(C_num, grid){
  
  .Call(C_cache_grid, C_num, grid$grid_data$x,  grid$grid_data$s_upper, grid$grid_data$p_a, grid$grid_data$s_upper_lower, grid$areas, grid$steps_number, grid$sampling_probabilities, grid$unif_scaler, grid$lt_properties, grid$rt_properties)
  
}

#' Grid builder
#'
#' @param lb scaler density lower bound
#' @param rb scaler density upper bound
#' @param mode vector density modes
#' @param f density function
#' @param h log transform of the density function
#' @param h_prime first derivative of h
#' @return message
#' @export
set_grid <- function(lb = -Inf, rb = Inf, mode, f = f_norm, h = h_norm, h_prime = h_prime_norm, to = c("set"), dist_name = "") {
  to <- match.arg(to)

  data_dir <- tools::R_user_dir("stors", "data")

  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

  grids_list_path <- file.path(data_dir, "grids_list.rds")

  if (!file.exists(grids_list_path)) {
    grids_list <- data.frame(distrubution_name = character(), number_of_steps = double(), sampling_efficiency = double())
  } else {
    grids_list <- readRDS(grids_list_path)

    stopifnot("Error: dist_name is already exist" = !(dist_name %in% grids_list$distrubution_name))
  }

  if (to == "set") {
    A <- seq(from = 0.005, to = 0.0002, length.out = 100)

    for (a in A) {
      opt_grid <- grid_builder(lb, rb, a, th = 0.6, mode, f, h, h_prime)

      if ((object.size(opt_grid) / 1024) > 128) break
    }


    func_to_text <- deparse(f)

    opt_grid$dens_func <- func_to_text

    new_grid_path <- file.path(data_dir, paste0(dist_name, ".rds"))

    saveRDS(opt_grid, file = new_grid_path)

    # grids_list = grids_list[nrow(grids_list) + 1,]=c(dist_name, opt_grid$steps_number, 1 / sum(opt_grid$areas) )
    grids_list <- rbind(grids_list, data.frame(distrubution_name = dist_name, number_of_steps = opt_grid$steps_number, sampling_efficiency = 1 / sum(opt_grid$areas)))

    print(grids_list)

    saveRDS(grids_list, grids_list_path)
  }
}


#' Title
#'
#' @return
#' @export
#'
#' @examples
grids_list_print <- function() {
  data_dir <- tools::R_user_dir("stors", "data")

  grids_list_path <- file.path(data_dir, "grids_list.rds")


  stopifnot(file.exists(grids_list_path))

  grids_list <- readRDS(grids_list_path)


  print(grids_list)
}




#' Title
#'
#' @param dist_name
#'
#' @return
#' @export
#'
#' @examples
read_grid <- function(dist_name) {
  data_dir <- tools::R_user_dir("stors", "data")

  stopifnot(`there are no grids in user_data directory` = dir.exists(data_dir))

  grid_dir <- file.path(data_dir, paste0(dist_name, ".rds"))

  stopifnot(`there is no grid for 'grid_dir' density` = file.exists(grid_dir))

  grid <- readRDS(grid_dir)

  grid$dens_func <- eval(parse(text = grid$dens_func))

  return(grid)
}
#' Grid init
#'
#' @param lb scaler density lower bound
#' @param rb scaler density upper bound
#' @param a scaler step area
#' @param th scaler pre_acceptance threshold
#' @param mode vector density modes
#' @param f density function
#' @param h log transform of the density function
#' @param h_prime first derivative of h
#' @return list including proposal distribution properties
#' @importFrom utils head
grid_builder <- function(lb = -Inf, rb = Inf, a, th, mode, f, h, h_prime) {
  mode_n <- length(mode)

  final_grid <- data.frame(
    x = c(), s_upper = c(), s_lower = c(), p_a = c(),
    s_upper_lower = c()
  )
  grids <- list()
  area <- c(0, 0, 0)
  g_len <- c()

  for (mode_i in (1:mode_n)) {
    grids[[mode_i]] <- find_steps(lb = -Inf, rb = Inf, a, th, mode[mode_i], mode_i, mode_n, f, h_prime, h)
  }


  area[1] <- grids[[1]]$l_tail_area



  if (mode_n > 1) {
    for (i in (1:(mode_n - 1))) {
      if (grids[[i]]$d$x[grids[[i]]$m + 1] > grids[[i + 1]]$d$x[1]) {
        if (grids[[i]]$d$s_upper[grids[[i]]$m] > grids[[i + 1]]$d$s_upper[1]) {
          grids[[i]]$d <- head(grids[[i]]$d, -1)

          grids[[i]]$d$s_upper[grids[[i]]$m] <- a / (grids[[i + 1]]$d$x[1] - grids[[i]]$d$x[grids[[i]]$m])

          grids[[i]]$d$p_a[grids[[i]]$m] <- 0

          grids[[i + 1]]$d$p_a[1] <- 0
        } else {
          grids[[i + 1]]$d$x[1] <- grids[[i]]$d$x[grids[[i]]$m + 1]

          grids[[i]]$d <- head(grids[[i]]$d, -1)


          grids[[i + 1]]$d$s_upper[1] <- a / (grids[[i + 1]]$d$x[2] - grids[[i + 1]]$d$x[1])

          grids[[i + 1]]$d$p_a[1] <- 0

          grids[[i]]$d$p_a[grids[[i]]$m] <- 0
        }
      } else {
        grids[[i]]$m <- grids[[i]]$m + 1

        grids[[i]]$d$s_upper[grids[[i]]$m] <- a / (grids[[i + 1]]$d$x[1] - grids[[i]]$d$x[grids[[i]]$m])

        grids[[i]]$d$p_a[grids[[i]]$m] <- 0
      }

      m1 <- grids[[i]]$m

      m2 <- grids[[i + 1]]$m




      area[2] <- area[2] + m1 * a

      g_len[length(g_len) + 1] <- m1



      if (i == (mode_n - 1)) {
        final_grid <- rbind(final_grid, grids[[i]]$d, grids[[i + 1]]$d)

        area[2] <- area[2] + m2 * a

        g_len[length(g_len) + 1] <- m2
      } else {
        final_grid <- rbind(final_grid, grids[[i]]$d)
      }
    }
  } else {
    final_grid <- grids[[1]]$d

    area[2] <- grids[[1]]$m * a

    g_len[length(g_len) + 1] <- grids[[1]]$m
  }

  area[3] <- grids[[mode_n]]$r_tail_area

  normalizing_con <- sum(area)

  area_cum_sum <- cumsum(area)

  sampling_probabilities <- (area_cum_sum / normalizing_con)[1:2]

  steps_number <- sum(g_len) # m

  unif_scaler <- normalizing_con / area[2]

  x1 <- final_grid$x[1]
  xm <- final_grid$x[steps_number + 1]

  lt_properties <- c(exp(h_upper(x1, lb, h_prime, h)), normalizing_con * h_prime(x1), h(x1), 1 / h_prime(x1), h_prime(x1))
  rt_properties <- c(normalizing_con, area_cum_sum[2], h_prime(xm) / f(xm), 1 / h_prime(xm), h_prime(xm), h(xm))

  invisible(list(grid_data = final_grid, areas = area, steps_number = steps_number, sampling_probabilities = sampling_probabilities, unif_scaler = unif_scaler, lt_properties = lt_properties, rt_properties = rt_properties))
}

#' Steps Builder
#'
#' @param lb scaler density lower bound
#' @param rb scaler density upper bound
#' @param a scaler step area
#' @param th scaler pre_acceptance threshold
#' @param mode vector density modes
#' @param mode_i mode index
#' @param mode_n total number of modes
#' @param f density function
#' @param h log transform of the density function
#' @param h_prime first derivative of h
#'
#' @return
#'
#' @examples
find_steps <- function(lb = -Inf, rb = Inf, a, th, mode, mode_i, mode_n, f, h_prime, h) {
  memory_res <- (max(500, ceiling(1 / a)) + 500) / 2

  x <- rep(NA, memory_res * 2 + 1)
  s_upper <- rep(NA, memory_res * 2 + 1)
  s_lower <- rep(NA, memory_res * 2 + 1)
  p_a <- rep(NA, memory_res * 2 + 1)
  s_upper_lower <- rep(NA, memory_res * 2 + 1)


  r <- 0
  l <- 0

  i <- memory_res + 1


  if (mode != rb) {
    x_c <- mode

    f_x <- f(x_c)

    while (TRUE) {
      x_next <- x_c + a / f_x

      if (x_next > rb || (mode_i == mode_n && f(x_next) / f_x < th)) {
        x[i + r] <- x_c
        if (rb == Inf) {
          r_tail_area <- (1 / h_prime(x[i + r])) * -f(x[i + r])
        } else {
          r_tail_area <- (1 / h_prime(x[i + r])) * (exp(h_upper(x[i + r], rb, h_prime, h)) - f(x[i + r]))
        }


        break
      }

      f_x_next <- f(x_next)

      if (f_x_next > f_x) {
        x[i + r] <- x_c
        s_upper[i + r] <- f_x
        r_tail_area <- 0
        break
      }


      x[i + r] <- x_c
      s_upper[i + r] <- f_x
      s_lower[i + r] <- f_x_next
      s_upper_lower[i + r] <- s_upper[i + r] / s_lower[i + r]
      p_a[i + r] <- s_lower[i + r] / s_upper[i + r]


      f_x <- f_x_next
      x_c <- x_next

      r <- r + 1
    }
  }



  if (mode != lb) {
    x_previous <- mode
    f_x_previous <- f(x_previous)


    while (TRUE) {
      x_c <- x_previous - a / f_x_previous



      if (x_c < lb || (mode_i == 1 && f(x_c) / f_x_previous < th)) {
        if (lb == Inf) {
          l_tail_area <- (1 / h_prime(x[i - l])) * f(x[i - l])
        } else {
          l_tail_area <- (1 / h_prime(x[i - l])) * (f(x[i - l]) - exp(h_upper(x[i - l], lb, h_prime, h)))
        }
        break
      }

      f_x <- f(x_c)

      if (f(x_c) > f_x_previous) {
        l_tail_area <- 0
        break
      }

      l <- l + 1
      x[i - l] <- x_c
      s_upper[i - l] <- f_x_previous
      s_lower[i - l] <- f_x
      s_upper_lower[i - l] <- s_upper[i - l] / s_lower[i - l]

      p_a[i - l] <- f_x / f_x_previous


      f_x_previous <- f_x
      x_previous <- x_c
    }
  }

  m <- l + r

  d <- data.frame(
    x = x, s_upper = s_upper, p_a = p_a,
    s_upper_lower = s_upper_lower
  )

  d <- subset(d, rowSums(is.na(d)) != ncol(d))

  return(list(d = d, m = m, l_tail_area = l_tail_area, r_tail_area = r_tail_area))
}





h_upper <- function(grid_point, val, h_prime, h) {
  h_prime(grid_point) * (val - grid_point) + h(grid_point)
}
