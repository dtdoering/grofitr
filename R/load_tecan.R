#' Load and clean data from BMG plate readers
#'
#' Loads in data from Tecan plate readers in a tidy (long) format.
#' Users can then directly plot the data or estimate growth parameters
#' with `grofitr()`.
#'
#' @param file Raw-text file output from Tecan plate reader containing growth data
#' @param get.barcode If TRUE, will search for a plate barcode in the "ID1"
#' field of the file header.
#'
#' @importFrom readr read_csv
#' @importFrom dplyr slice mutate rename_at vars select %>% funs
#' @importFrom tidyr gather
#'
#' @export

# Add option to extract barcode? Would store in "plate" column.

load_tecan <- function(file, get.barcode = FALSE) {
  x <- file %>%
    readLines() %>%
    grepl("^0s", .) %>%
    which() %>%
    `-`(2) %>%
    read_csv(file, skip = .)
  y <- x %>% rename(time = X1, temp = X2) %>%
    dplyr::slice(-min(which(grepl("^Date", time))):-n()) %>%
    tidyr::gather(well, OD, -temp, -time) %>%
    dplyr::mutate(time = time %>% str_replace("s", "") %>% as.numeric() %>% magrittr::divide_by(3600))
  if (get.barcode == TRUE) {
    y <-  dplyr::mutate(.data = y, plate = get_barcode_BMG(file))
  } else {
    y <-  dplyr::mutate(.data = y, plate = basename(file) %>%
                          strsplit("\\.") %>% `[[`(1) %>% `[`(1))
  }
}
