#'
#' @importFrom magrittr %>%
#' @importFrom stringr str_extract

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
