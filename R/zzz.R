#' @useDynLib stors, .registration = TRUE, .fixes = "C_"
# Globals:


built_in_proposals <- list(
  srnorm = list(
    c_num = 1,
    tails_method = "ARS",
    scalable = TRUE,
    std_params = list(mean = 0, sd = 1),
    transform_params = function(par) {
      par$sd <- 1 / par$sd
      return(par)
    },
    create_f = function(mu, sd) {
      fun_txt <- paste0(
        "function(x){((1.0 / (",
        sd, ")*exp(-0.5 * ((x -"
        , mu, ")/", sd,
        ")*((x - ", mu,
        ") / ", sd, "))))}")
      return(eval(parse(text = fun_txt)))
    },
    set_modes = function(mu = 0) {
      return(mu)
    },
    lb = -Inf,
    rb = Inf
  ), srlaplace = list(
    c_num = 3,
    tails_method = "IT",
    scalable = TRUE,
    std_params = list(mu = 0, b = 1),
    transform_params = function(par) {
      return(par)
    },
    create_f = function(mu, b) {
      fun_txt <- paste0("function(x) { (1 / (2 *", b
                        , ")) * exp(-abs(x -", mu,
                        ") /", b, ")}")
      return(eval(parse(text = fun_txt)))
    },
    create_cdf = function(mu, b) {
      function(x) {
        if (x <= mu) {
          0.5 * exp((x - mu) / b)
        } else {
          1 - 0.5 * exp(-(x - mu) / b)
        }
      }
    },
    set_modes = function(mu = 0) {
      mu
    },
    lb = -Inf,
    rb = Inf
  ), srexp = list(
    c_num = 5,
    tails_method = "IT",
    scalable = TRUE,
    std_params = list(rate = 1),
    transform_params = function(par) {
      return(par)
    },
    create_f = function(rate) {
      fun_txt <- paste0("function(x) { fx <-", rate, " * exp(-", rate, " * x)
                         fx <- ifelse(is.nan(fx), 0, fx)
                         return(fx) }")
      return(eval(parse(text = fun_txt)))
    },
    create_cdf = function(rate) {
      function(x) {
        return(1 - exp(- rate * x))
      }
    },
    lb = 0,
    rb = Inf
  ), srchisq = list(
    c_num = 7,
    tails_method = "ARS",
    scalable = FALSE,
    transform_params = function(par) {
      return(par)
    },
    create_f = function(df) {
      fun_txt <- paste0(
        "function(x) {
       fx <- (1 / (2 ^ (",
        df, " / 2) * gamma(",
        df, " / 2))) * x ^ (",
        df, "/ 2 - 1) * exp(-x / 2)
        fx <- ifelse(is.nan(fx), 0, fx)
       return(fx)}"
      )
      return(eval(parse(text = fun_txt)))
    },
    set_modes = function(df = 1) {
      max(df - 2, 0)
    },
    lb = 0,
    rb = Inf
  ), srgamma = list(
    c_num = 9,
    std_params = list(shape = 1, scale = 1),
    tails_method = "ARS",
    scalable = FALSE,
    transform_params = function(par) {
      return(par)
    },
    create_f = function(shape = 1, scale = 1) {
      rate <- 1 / scale
      fun_txt <- paste0(
        "function(x) {
        fx <- ifelse(x < 0, 0, (1 / (gamma(", shape, ") * ", scale, "^", shape, ") * x ^ (", shape, " - 1) * exp(-x /", scale, ")))
        fx <- ifelse(is.nan(fx), 0, fx)
        return(fx)
      }"
      )
      return(eval(parse(text = fun_txt)))
    },
    set_modes = function(shape, scale) {
      if (shape < 1) 0 else (shape - 1) * scale
    },
    lb = 0,
    rb = Inf
  ),  srbeta = list(
    c_num = 11,
    tails_method = "ARS",
    scalable = FALSE,
    std_params = list(shape1 = 2, shape2 = 2),
    transform_params = function(par) {
      return(par)
    },
    create_f = function(shape1, shape2) {
      fun_txt <- paste0(
        "function(x){
        fx <- (x^(", shape1, "-1) * (1 - x)^(", shape2, " - 1)) / beta(", shape1, ",", shape2, ")
        return(fx)
        }")
      return(eval(parse(text = fun_txt)))
    },
    set_modes = function(shape1, shape2) {
      return((shape1 - 1) / (shape1 + shape2 - 2))
    },
    lb = 0,
    rb = 1
  ), srpareto = list(
    c_num = 13,
    tails_method = "IT",
    scalable = FALSE,
    std_params = list(scale = 1, shape = 1),
    transform_params = function(par) {
      return(par)
    },
    create_f = function(scale, shape) {
      fun_txt <- paste0("function(x) { fx <- ifelse(x < ", scale, ", 0, (", shape * scale ^ shape, ") / x^( ",
                        shape, " + 1 ))", "
                        return(fx) }")
      return(eval(parse(text = fun_txt)))
    },
    create_cdf = function(scale, shape) {
      function(x) {
        return(1 - (scale / x) ^ shape)
      }
    },
    set_modes = function(scale) {
      return(scale)
    },
    lb = NULL,
    rb = Inf
  )
)

stors_env <- new.env(parent = emptyenv())


.onLoad <- function(lib, pkg) {

  data_dir <- tools::R_user_dir("stors", "data")

  if (!dir.exists(data_dir))
    dir.create(data_dir, recursive = TRUE)

  if (!file.exists(file.path(data_dir, "version"))) {
    # No versioning file so make sure directory is empty and create new version file
    unlink(list.files(tools::R_user_dir("stors", "data"), full.names = TRUE),
           recursive = TRUE, force = TRUE)
    cat(paste(utils::packageVersion("stors"), sep = "."),
        file = file.path(data_dir, "version"))
  } else {
    vers <- suppressWarnings(readLines(file.path(data_dir, "version"), n = 1))
    if (!identical(paste(utils::packageVersion("stors"), sep = "."),
                   vers)) {
      warning("Package version updated, old grids being archived.")
      if (file.exists(file.path(data_dir, paste0("builtin_grids_", vers)))) {
        unlink(file.path(data_dir, c(paste0("builtin_grids_", vers), paste0("user_proposals_", vers))), recursive = TRUE, force = TRUE)
      }
      file.rename(file.path(data_dir, "builtin_grids"),
                  file.path(data_dir, paste0("builtin_grids_", vers)))
      file.rename(file.path(data_dir, "user_proposals"),
                  file.path(data_dir, paste0("user_proposals_", vers)))
      cat(paste(utils::packageVersion("stors"), sep = "."),
          file = file.path(data_dir, "version"))
    }
  }

  builtin_grids_dir <- file.path(data_dir, "builtin_grids")

  user_proposals_dir <- file.path(data_dir, "user_proposals")


  if (!dir.exists(builtin_grids_dir))
    dir.create(builtin_grids_dir, recursive = TRUE)


  if (!dir.exists(user_proposals_dir))
    dir.create(user_proposals_dir, recursive = TRUE)

  user_cnum_counter <- 101
  user_session_cached_grid_locks <- data.frame(lock = character(), cnum = numeric())

  assign("builtin_grids_dir", builtin_grids_dir, envir = stors_env)
  assign("user_proposals_dir", user_proposals_dir, envir = stors_env)
  assign("user_cnum_counter", user_cnum_counter, envir = stors_env)
  assign("user_session_cached_grid_locks", user_session_cached_grid_locks, envir = stors_env)


  builtin_grids <- list.files(builtin_grids_dir)

  if (length(builtin_grids) == 0) {


    for (name in names(built_in_proposals)) {

      fun_name <- paste0(name, "_optimize")
      do.call(fun_name, list(steps = 4091))

    }

  } else {

    for (proposal_name in builtin_grids) {

      proposal_path <- file.path(builtin_grids_dir, proposal_name)
      grid <- readRDS(proposal_path)

      if (is_valid_proposal(grid))
        cache_proposal_c(grid$cnum, grid)



    }

  }

}

.onUnload <- function(...) {
  .Call(C_free_cache)
}
