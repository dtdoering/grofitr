#' Extract plate barcodes from BMG plate reader data files
#'
#' Looks for a barcode present in the "ID1" field of plate reader data files
#' output by BMG Labtech instruments, such as FLUOstar, SPECTROstar, etc.
#'
#' @param file File path to the plate whose barcode you wish to extract
#' @param trim If TRUE, trim leading 0s off of the extracted barcode (default = FALSE).
#' @importFrom magrittr %>%
#' @importFrom stringr str_extract
#'
#' @export

get_barcode_BMG <- function(file, trim = FALSE) {
  if (trim == TRUE) {
    file %>%
      readLines() %>%
      grep("ID1: ", ., value = T) %>%
      str_extract("(?<=ID1: ).*(?= \\d{1,2}/\\d{1,2}/\\d{2,4})") %>%
      substr(regexpr("[^0]", .), nchar(.))
  } else if (trim == FALSE) {
    file %>%
      readLines() %>%
      grep("ID1: ", ., value = T) %>%
      str_extract("(?<=ID1: ).*(?= \\d{1,2}/\\d{1,2}/\\d{2,4})")
  }
}
