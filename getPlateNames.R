options(stringsAsFactors = FALSE, warn = 0)

library("magrittr")

getPlateNames <- function(path) {
  # Create file names for the saved plate reader data
  DataFiles <-
  list.files(path) %>% grep("(-E\\d-P\\d|TRno\\d*)\\.(CSV|csv)", ., value = T)

  # Create a vector of object names using barcodes found in csv file
  PlateNames <- character() # Will become 4-digit barcodes
  for (i in seq_along(DataFiles)) {
      PlateNames[i] <- grep("ID1: ",
                            path %>% paste(DataFiles[i], sep = "") %>% readLines(warn = F),
                            value = T) %>%
                       substr(6,12) %>%
                       substr(regexpr("[^0]", .),
                              nchar(.))
  }
  return(PlateNames)
}
