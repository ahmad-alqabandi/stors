#' @useDynLib stors, .registration = TRUE, .fixes = "C_"
# Globals:


built_in_proposals <- list(
  srnorm = list(
    name = "srnorm",
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
    lower = -Inf,
    upper = Inf
  ), srlaplace = list(
    name = "srlaplace",
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
    lower = -Inf,
    upper = Inf
  ), srexp = list(
    name = "srexp",
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
    lower = 0,
    upper = Inf
  ), srchisq = list(
    name = "srchisq",
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
    lower = 0,
    upper = Inf
  ), srgamma = list(
    name = "srgamma",
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
    lower = 0,
    upper = Inf
  ),  srbeta = list(
    name = "srbeta",
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
    lower = 0,
    upper = 1
  ), srpareto = list(
    name = "srpareto",
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
    lower = NULL,
    upper = Inf
  )
)

stors_env <- new.env(parent = emptyenv())


.onLoad <- function(lib, pkg) {

  data_dir <- tools::R_user_dir("stors", "data")

  if (!dir.exists(data_dir))
    dir.create(data_dir, recursive = TRUE)

  # if (!file.exists(file.path(data_dir, "version"))) {
  #   # No versioning file so make sure directory is empty and create new version file
  #   unlink(list.files(tools::R_user_dir("stors", "data"), full.names = TRUE),
  #          recursive = TRUE, force = TRUE)
  #   fd <- file(file.path(data_dir, "version"), open = "wt")
  #   write(as.character(utils::packageVersion("stors")), file = fd)
  #   close(fd)
  # } else {
  #   vers <- suppressWarnings(readLines(file.path(data_dir, "version"), n = 1))
  #   if (!identical(as.character(utils::packageVersion("stors")),
  #                  vers)) {
  #     warning("Package version updated, old proposals being archived.")
  #     if (file.exists(file.path(data_dir, paste0("builtin_proposals_", vers)))) {
  #       unlink(file.path(data_dir, c(paste0("builtin_proposals_", vers), paste0("user_proposals_", vers))), recursive = TRUE, force = TRUE)
  #     }
  #     file.rename(file.path(data_dir, "builtin_proposals"),
  #                 file.path(data_dir, paste0("builtin_proposals_", vers)))
  #     file.rename(file.path(data_dir, "user_proposals"),
  #                 file.path(data_dir, paste0("user_proposals_", vers)))
  #     fd <- file(file.path(data_dir, "version"), open = "wt")
  #     write(as.character(utils::packageVersion("stors")), file = fd)
  #     close(fd)
  #   }
  # }

  builtin_proposals_dir <- file.path(data_dir, "builtin_proposals")

  user_proposals_dir <- file.path(data_dir, "user_proposals")


  if (!dir.exists(builtin_proposals_dir))
    dir.create(builtin_proposals_dir, recursive = TRUE)


  if (!dir.exists(user_proposals_dir))
    dir.create(user_proposals_dir, recursive = TRUE)

  user_cnum_counter <- 100
  user_session_cached_proposals_locks <- data.frame(lock = character(), cnum = numeric())

  assign("builtin_proposals_dir", builtin_proposals_dir, envir = stors_env)
  assign("user_proposals_dir", user_proposals_dir, envir = stors_env)
  assign("user_cnum_counter", user_cnum_counter, envir = stors_env)
  assign("user_session_cached_proposals_locks", user_session_cached_proposals_locks, envir = stors_env)


  builtin_proposals_files <- list.files(builtin_proposals_dir)
  existing_builtin_proposals_number <- as.numeric(sub("\\.rds$", "", builtin_proposals_files))
  number_of_proposals <- 2 * length(built_in_proposals)

  for (i in seq_len(number_of_proposals)) {

    proposal_path <- if (i %in% existing_builtin_proposals_number) {
      file.path(builtin_proposals_dir, paste0(i, ".rds"))
    } else {
      system.file(paste0("builtin_proposals/", i, ".rds"), package = pkg)
    }

    if (file.exists(proposal_path)) {

      proposal <- readRDS(proposal_path)

      if (is_valid_proposal(proposal)) {
        cache_proposal_c(proposal$cnum, proposal)
      }
    }

  }

}

.onUnload <- function(...) {
  .Call(C_free_cache)
}
