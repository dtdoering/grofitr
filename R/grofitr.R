#'
#' @import magrittr
#' @import readr
#' @import dplyr
#' @import tidyr
#'

grofitr <- function(plate, ...) {
  timepoints <- plate %>% select(time) %>% unique() %>% unlist()
  n <- length(timepoints)
  times <- timepoints %>% rep(96) %>% matrix(c(n, 96)) %>% t()
  dots <- list(...)

  plate <- plate %>%
    dplyr::mutate(OD = as.numeric(OD)) %>%
    tidyr::spread(time, OD) %>%
    dplyr::select(-Content) %>%
    dplyr::select(3,1,2,4:length(.))

  grofit::grofit(times, plate, ec50 = F,
                 control = grofit::grofit.control(neg.nan.act = F, interactive = F,
                                                  suppress.messages = T, ...)) %$% gcFit %$% gcTable %>% rename(`Well Row` = Well.Row, `Well Col` = Well.Col)
}
