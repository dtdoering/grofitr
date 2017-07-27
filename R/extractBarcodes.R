#' @export

options(stringsAsFactors = FALSE, warn = 0)

extract.barcodes <- function(path) {
  # Create file names for the saved plate reader data
  data.files <-
  list.files(path) %>% grep("(-E\\d-P\\d|TRno\\d*)\\.(CSV|csv)", ., value = T)

  # Create a vector of object names using barcodes found in csv file
  barcodes <- character() # Will become 4-digit barcodes
  for (i in seq_along(data.files)) {
      barcodes[i] <- grep("ID1: ",
                          path %>% paste(data.files[i], sep = "") %>%
                            readLines(warn = F),
                          value = T) %>%
                       substr(6, 12) %>%
                       substr(regexpr("[^0]", .),
                              nchar(.))
  }
  return(barcodes)
}
