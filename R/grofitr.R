#' @importFrom magrittr %>% %$% extract2
#' @importFrom dplyr mutate select rename
#' @importFrom tidyr spread
#' @import grofit
#'

grofitr <- function(plate, ...) {
  timepoints <- plate %>% select(time) %>% unique() %>% unlist()
  n <- length(timepoints)
  times <- timepoints %>% rep(96) %>% matrix(c(n, 96)) %>% t()

  plate <- plate %>%
    mutate(OD = as.numeric(OD)) %>%
    spread(time, OD) %>%
    select(-Content) %>%
    select(3,1,2,4:length(.))

  grofit(times,
         plate,
         ec50 = F,
         control = grofit::grofit.control(neg.nan.act = F,
                                          interactive = F,
                                          suppress.messages = T,
                                          ...)) %>%
         extract2("gcFit") %>%
         extract2("gcTable") %>%
         rename(`Well Row` = Well.Row,
                `Well Col` = Well.Col)
}
