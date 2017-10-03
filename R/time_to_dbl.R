#' Convert time to decimal
#'
#' Converts time in the form "HH h MM min" to hours in decimal form
#'
#' @param t A character vector of times in the form "HH h\[ MM min\]" or "X.Y"
#'
#' @return A numeric vector of times in hour-decimal form
#'
#' @examples
#' time_to_dbl("12 h 1 min")
#' t <- c("0 h", "0 h 45 min", "1 h 30 min")
#' time_to_dbl(t)
#'
#' @author Drew T. Doering, \email{dtdoering@@wisc.edu}
#'
#' @export
#' @importFrom magrittr %>% divide_by add
#' @importFrom stringr str_replace_all str_split

time_to_dbl <- function(t) {
  if ( grepl("h", t[1]) ) { # in "X h Y min" format
    t %>%
      str_replace_all("(\\d+) h( (\\d+) min)?", "\\1.0\\3") %>%
      str_split(string = ., pattern = "\\.") %>%
      sapply(FUN = function(x) {
        x %>% as.numeric() %>%
          `[`(2) %>%
          divide_by(60) %>%
          add(as.numeric(x[1]))
      })
  } else if ( any(grepl('\\.', t)) ) { # already in decimal format
    t %>% as.numeric()
  }
}
