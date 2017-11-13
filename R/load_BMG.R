#' Load and clean data from BMG plate readers
#'
#' Loads in data from BMG Labtech plate readers such as FLUOstar, CLARIOstar,
#' etc. in a tidy (long) format. Users can then directly plot the data or
#' estimate growth parameters with `grofitr()`.
#'
#' @param file Raw-text file output from BMG plate reader containing growth data
#' @param get.barcode If TRUE, will search for a plate barcode in the "ID1"
#' field of the file header.
#'
#' @importFrom magrittr %>% subtract
#' @importFrom readr read_csv
#' @importFrom dplyr slice mutate rename_at vars select
#' @importFrom tidyr gather

# Add option to extract barcode? Would store in "plate" column.

load_BMG <- function(file, get.barcode = FALSE) {
  x <- file %>%
    readLines() %>%
    grepl("Well Row", .) %>%
    which() %>%
    subtract(1) %>%
    read_csv(file, skip = .)
  y <- x %>% rename_at(vars(-c(1:3)),
                  funs(x %>% slice(1) %>% select(-c(1:3)))) %>%
    slice(-1) %>%
    gather(time, OD, -`Well Row`, -`Well Col`, -Content) %>%
    mutate(time = time_to_dbl(time), OD = as.numeric(OD))
  if (get.barcode == TRUE) {
    y %>% mutate(Barcode = get_barcode_BMG(file))
  } else {
    y %>% mutate(plate = basename(file) %>%
                        strsplit("\\.") %>% `[[`(1) %>% `[`(1))
  }
}
