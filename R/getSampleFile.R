getSampleFile <- function(filelist, animalId, cellName) {
  files <- grep(pattern = glue::glue("{animalId}/"), x = filelist, value = TRUE)
  files <- grep(pattern = glue::glue("{cellName}/"), x = files, value = TRUE)

  stopifnot(length(files) == 1)

  files
}
