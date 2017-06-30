library(readxl)

dataset.paths <- list.files("data-raw", "*.xls")
dats <- lapply(dataset.paths, function(path) {
  dat <- as.data.frame(read_excel(file.path("data-raw", path)))
  colnames(dat) <- tolower(colnames(dat))
  # Log-transform
  dat <- log(dat)
  dat$source <- path
  dat
})
flowCytometry <- Reduce(rbind, dats)

use_data(flowCytometry)
