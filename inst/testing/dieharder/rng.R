
set.seed(212)
ssize = 100000000

name = "runif_data"

runif_data <- as.integer(runif(ssize)*2^32-2147483648)

# Get the current working directory
current_dir <- paste0(getwd(),"/inst/testing/dieharder")

# Specify the file path in the current directory
file_path <- file.path(current_dir, paste0(name, ".bin"))

# Save the vector as a binary file
writeBin(runif_data, file_path, size = 4)

set.seed(212)
name = "rnorm_data"

rnorm_data <- rnorm(ssize)
rnorm_data <- (rnorm_data - min(rnorm_data))/(max(rnorm_data)-min(rnorm_data))
rnorm_data <- as.integer(rnorm_data*(2^32-2)-2147483647)

# Get the current working directory
current_dir <- paste0(getwd(),"/inst/testing/dieharder")

# Specify the file path in the current directory
file_path <- file.path(current_dir, paste0(name, ".bin"))

# Save the vector as a binary file
writeBin(runif_data, file_path, size = 4)







writeLines(c("type: d", "count: 1000000", "numbit: 32", runif(1000000)*2^32), paste0(getwd(),"/inst/testing/dieharder/runif_data.input"))
writeLines(c("type: d", "count: 1000000", "numbit: 32", rep(2024, 1000000)), paste0(getwd(),"/inst/testing/dieharder/const_data.input"))
dat <- rnorm(1000000)
dat <- (dat-min(dat))/(max(dat)-min(dat))*2^32
writeLines(c("type: d", "count: 1000000", "numbit: 32", round(dat)), paste0(getwd(),"/inst/testing/dieharder/rnorm_data.input"))


name = "srnorm_data"


data <- srnorm(ssize)

# Get the current working directory
current_dir <- paste0(getwd(),"/inst/testing/dieharder")

# Specify the file path in the current directory
file_path <- file.path(current_dir, paste0(name, ".bin"))

# Save the vector as a binary file
writeBin(data, file_path)



name = "rnorm_data"

data <- rnorm(ssize)

# Get the current working directory
current_dir <- paste0(getwd(),"/inst/testing/dieharder")

# Specify the file path in the current directory
file_path <- file.path(current_dir, paste0(name, ".bin"))

# Save the vector as a binary file
writeBin(data, file_path)




name = "zrnorm_data"

data <- zrnormR(ssize)

# Get the current working directory
current_dir <- paste0(getwd(),"/inst/testing/dieharder")

# Specify the file path in the current directory
file_path <- file.path(current_dir, paste0(name, ".bin"))

# Save the vector as a binary file
writeBin(data, file_path)

