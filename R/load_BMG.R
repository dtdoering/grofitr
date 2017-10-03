#'
#' @import magrittr
#' @import readr
#' @import dplyr
#' @import tidyr

# Add option to extract barcode? Would store in "plate" column.

load_BMG <- function(file, get.barcode = FALSE) {
  x <- file %>%
    readLines() %>%
    grepl("Well Row", .) %>%
    which() %>%
    subtract(1) %>%
    read_csv(file, skip = .)# %>%
  y <- x %>% rename_at(vars(-c(1:3)),
                  funs(x %>% slice(1) %>% select(-c(1:3)))) %>%
    dplyr::slice(-1) %>%
    tidyr::gather(time, OD, -`Well Row`, -`Well Col`, -Content) %>%
    dplyr::mutate(time = time_to_dbl(time), OD = as.numeric(OD))
  if (get.barcode == TRUE) {
    y %>%  dplyr::mutate(Barcode = get_barcode_BMG(file))
  } else {
    y %>%  dplyr::mutate(plate = basename(file) %>%
                        strsplit("\\.") %>% `[[`(1) %>% `[`(1))
  }
}
