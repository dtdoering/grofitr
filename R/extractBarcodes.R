#' INPUT = full file path
#' OUTPUT = Barcode extracted from ID1 field, if there is one there
#'
#' @import magrittr
#' @import readr
#' @import dplyr
#' @import tidyr


extract_barcode <- function(file) {

  # Create file names for the saved plate reader data
  data.files <- list.files(file) %>%
    grep("(-E\\d-P\\d|TRno\\d*)\\.(CSV|csv)", ., value = T)

  # Create a vector of object names using barcodes found in csv file
  barcodes <- character(length = 1) # Will become 4-digit barcodes
  for (i in seq_along(data.files)) {
      barcodes[i] <- grep("ID1: ",
                          file %>% paste(data.files[i], sep = "") %>%
                            readLines(warn = F),
                          value = T) %>%
                       substr(6, 12) %>%
                       substr(regexpr("[^0]", .),
                              nchar(.))
  }
  return(barcodes)
}
