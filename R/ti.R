#' @title Theil Index (TI) and Decomposition

#' @description Calculate the theil index, and decompose it between groups.

#' @param y Outcome varialbe.
#' @param s Variable based to rank y.
#' @param group1 First level of cases.
#' @param group2 Second level of cases.

#' @author Shanquan CHEN \email{shanquan0301@gmial.com}

#' @keywords Theil Index
#' @examples
##' library(dineq)
##' data(mex_inc_2016)
##' mex_inc_2016 %$% ti(y = age,
##'                  s = income,
##'                  group1 = hh_number,
##'                  group2 = hh_structure)



#' @import dplyr magrittr

#' @export ti
#'
ti <- function(y = y, s = s, group1 = group1, group2 = group2){
  res <- data.frame("ti" = numeric(1), "ti_between" = numeric(1), "ti_within" = numeric(1),
                    "ti_between_rate" = numeric(1), "ti_within_rate" = numeric(1))
  y[y == 0] <- ifelse(length(unique(y)) == 2, 0.1, 0.1)
  s[s == 0] <- 0.01
  mdat1 <- data.frame(y, s, group1, group2)
  y_total <- mdat1 %$% sum(y, na.rm = TRUE)
  s_total <- mdat1 %$% sum(s, na.rm = TRUE)
  mdat2 <- mdat1 %>% group_by(group2) %>%
    summarise(y_sub_total = sum(y, na.rm = TRUE),
              s_sub_total = sum(s, na.rm = TRUE),
              y_b = log10((y_sub_total/y_total) / (s_sub_total/s_total)) * (y_sub_total/y_total),
              y_w = sum(log10((y/y_sub_total) / (s/s_sub_total)) * (y/y_total), na.rm = TRUE))
  theil_bet <- sum(mdat2$y_b, na.rm = TRUE)
  theil_win <- sum(mdat2$y_w, na.rm = TRUE)
  res$ti <- theil_bet + theil_win
  res$ti_between <- theil_bet
  res$ti_within <- theil_win
  res$ti_between_rate <- res$ti_between * 100 / res$ti
  res$ti_within_rate <- res$ti_within * 100 / res$ti
  return(res)
}
