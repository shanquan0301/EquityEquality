#' @title Horizontal Index (HI) and Decomposition

#' @description Calculate the horizontal index, and decompose it among factors.

#' @param inc Variable based to rank y.
#' @param y is the dependent variable.
#' @param x is the data frame of independent variables.
#' @param fir.st if whether fit for first step. Default is TRUE.
#' @param nam_adj variables needed to be transfered into dummy variavles. Variables list in nam_adj must be factors.
#' @param need_fa variables used to reflect the health need.


#' @author Shanquan CHEN \email{shanquan0301@gmial.com}

#' @keywords Horizontal Index
#' @examples
##' res <- dat_2011_rural_imp %$%
##'             fun_hi_sd2(inc = hous_inc,
##'             y = inp_yes,
##'             x = data.frame(health_professionals.10.000., te,
##'                          edu, live_alone, nation, smoke_year, drink_year, phy_act, liv_con, mental_num, fam_n_inc,
##'                          hous_inc, ncmi, pens, comi,
##'                          age_dec, gender, self_stat, adl_heal, pain_num, chro_num),
##'             nam_adj = c("edu", "self_stat"),
##'             need_fa = c("age_dec", "gender", "self_stat", "adl_heal", "pain_num", "chro_num"))

#' @import dplyr magrittr stringr stats mfx MASS car
#' @importFrom tidyr spread
#' @importFrom nnet class.ind

#' @export fun_hi_sd2
#'
#----------------------------------------------
fun_hi_sd2 <- function(inc, y, x, fir.st = TRUE, nam_adj = nam_adj, need_fa = need_fa){
  rank.inc.kna <- rank(inc, na.last = "keep")
  rank.inc.rna <- rank(inc, na.last = NA)
  R.inc <- (rank.inc.kna - 0.5) / length(rank.inc.rna)
  var.R <- var(R.inc, na.rm = TRUE)
  res1 <- data.frame("mean" = numeric(0),"sd" = numeric(0),"x2.5 %" = numeric(0), "x97.5 %" = numeric(0))
  res2 <- data.frame("Estimate" = numeric(0), "Pr(>|t|)" = numeric(0), "ci2.5 %" = numeric(0), "ci97.5 %" = numeric(0))
  res_sep <- data.frame("mean" = numeric(0),"sd" = numeric(0),"x2.5 %" = numeric(0), "x97.5 %" = numeric(0),
                        "Estimate" = numeric(0), "Pr(>|t|)" = numeric(0), "ci2.5 %" = numeric(0), "ci97.5 %" = numeric(0))
  #---------
  #browser()
  for(i in 1: length(x)){
    if(sum(!x[, i] == 0) == 0){
      x[, i] <- 0.001
    }
  }
  #----------
  #---------------------------------
  #dummy the factor in x
  x_dam <- fun_x_dam2(data = x, nam_adj = nam_adj)
  #------------------------------
  #calculate the mean of non need factor
  x_mean <- x_dam
  for (i in names(x_mean)){
    if (sum(str_detect(i, need_fa)) == 0)
      x_mean[, i] <- mean(x_mean[, i])
  }
  #-------------------------------
  # sub the non need factor in x_dam
  key <- NA
  for (i in 1:length(names(x_dam))){
    if (sum(str_detect(names(x_dam)[i], need_fa)) > 0)
      key[i] <- FALSE
    else
      key[i] <- TRUE
  }
  x_non_need <- x_dam[, key]
  #-------------------------------------
  #browser()
  if (fir.st == TRUE) reg <- glm(y ~ ., data = x, family = binomial(link = "logit")) else reg <- glm.nb(y ~ ., data = x)
  y_mat  <- exp(predict(reg, new.data = x_mean)) / (1 + exp(predict(reg, new.data = x_mean)))
  y_is <- y - y_mat + mean(y_mat)
  x_mid <- cbind(y_is, x_non_need)
  for(i in 1:length(x_mid)){
    x.handle <- x_mid[, i]
    if (is.factor(x.handle) == TRUE) x.handle <- as.numeric(x.handle)
    reg2 <-  lm((2 * var.R / mean(x.handle, na.rm = TRUE)) * x.handle ~ R.inc)
    res1[i,1] <- mean(x.handle, na.rm = TRUE)
    res1[i,2] <- sd(x.handle, na.rm = TRUE)
    res1[i,3:4] <- quantile(x.handle, probs = c(0.25, 0.75), na.rm = TRUE)
    res2 <- rbind(res2,cbind(t(summary(reg2)$coefficients[2, c(1, 4)]), t(confint(reg2)[2, ])))
  }
  res_sep <- cbind(res1, res2)
  #------------------
  for(i in 1:length(x_mid)){
    if (res_sep[i,6] <= 0.001) res_sep$sig[i] <- "***"
    else
      if(res_sep[i,6] <= 0.01) res_sep$sig[i] <- "**"
      else
        if (res_sep[i,6] <= 0.05) res_sep$sig[i] <- "*"
        else
          if (res_sep[i,6] <= 0.1) res_sep$sig[i] <- "."
          else res_sep$sig[i] <- " "
  }
  #-------------------
  row.names(res_sep) <- c("y_HI", names(x_non_need))
  names(res_sep) <- c("mean", "sd", "x2.5 %", "x97.5 %", "Estimate", "Pr(>|t|)", "ci2.5 %", "ci97.5 %", "sig")
  #print(res_sep)
  # this part is for res_com
  res_com <- data.frame("mean(sd)" = character("0"), "CI" = character("0"))
  res_com$mean.sd. <- as.character(res_com$mean.sd.)
  res_com$CI <- as.character(res_com$CI)
  for (i in 1:length(x_mid)){
    res_com[i, 1] <- str_c(round(res_sep[i, 1], 4), "(", round(res_sep[i, 2], 4), ")")
    res_com[i, 2] <- str_c(round(res_sep[i, 5], 4), "(", round(res_sep[i, 7], 4), ", ", round(res_sep[i, 8], 4), ")")
  }
  row.names(res_com) <- row.names(res_sep)
  res_sep <- format(res_sep, justify = "left", digits = 4)
  #this part is for regression
  reg3 <- lm(y_is ~ ., data = x_non_need)
  #this part is for decomposition
  decom <- as.data.frame(summary(reg3)$coefficients[-1, ])
  decom$contri.coef <- decom$Estimate * as.numeric(res_sep$mean[-1]) / mean(y_is, na.rm = TRUE)
  decom$contri <- decom$contri.coef * as.numeric(res_sep$Estimate[-1])
  decom <- rbind(decom, c(0, 0, 0, NA, 0, 0))
  decom$contri[length(x_non_need) + 1] <- as.numeric(res_sep$Estimate[1]) - sum(decom$contri, na.rm = TRUE)
  decom$contri.rate <- decom$contri * 100 / abs(as.numeric(res_sep$Estimate[1])) # this rate is based on HI
  row.names(decom)[length(x_non_need) + 1] <- "residual"
  #-------------
  decom$sig <- car::recode(decom$`Pr(>|t|)`, "0:0.001 = '***'; 0.001:0.01 = '**'; 0.01:0.05 = '*'; 0.05:0.1 = '.'; else = ' '")
  decom <- format(decom, justify = "left", digits = 4)
  res <- list(res_sep, res_com, decom)
  names(res) <- c("res_sep", "res_com", "decom")
  return(res)
}
#------------------------------------------
#hi = CI_m - CI_n, only decom the CI_m
fun_hi_sd2.1 <- function(inc, y, x, fir.st = TRUE, nam_adj = nam_adj, need_fa = need_fa){
  rank.inc.kna <- rank(inc, na.last = "keep")
  rank.inc.rna <- rank(inc, na.last = NA)
  R.inc <- (rank.inc.kna - 0.5) / length(rank.inc.rna)
  var.R <- var(R.inc, na.rm = TRUE)
  res1 <- data.frame("mean" = numeric(0),"sd" = numeric(0),"x2.5 %" = numeric(0), "x97.5 %" = numeric(0))
  res2 <- data.frame("Estimate" = numeric(0), "Pr(>|t|)" = numeric(0), "ci2.5 %" = numeric(0), "ci97.5 %" = numeric(0))
  res_sep <- data.frame("mean" = numeric(0),"sd" = numeric(0),"x2.5 %" = numeric(0), "x97.5 %" = numeric(0),
                        "Estimate" = numeric(0), "Pr(>|t|)" = numeric(0), "ci2.5 %" = numeric(0), "ci97.5 %" = numeric(0))
  #---------
  for(i in 1: length(x)){
    if(sum(!x[, i] == 0) == 0)
      x[, i] <- 0.001
  }
  #----------
  #---------------------------------
  #dummy the factor in x
  x_dam <- fun_x_dam2(data = x, nam_adj = nam_adj)
  #------------------------------
  #calculate the mean of non need factor
  x_mean <- x_dam
  for (i in names(x_mean)){
    if (sum(str_detect(i, need_fa)) == 0)
      x_mean[, i] <- mean(x_mean[, i])
  }
  #-------------------------------
  # sub the non need factor in x_dam
  key <- NA
  for (i in 1:length(names(x_dam))){
    if (sum(str_detect(names(x_dam)[i], need_fa)) > 0)
      key[i] <- FALSE
    else
      key[i] <- TRUE
  }
  x_non_need <- x_dam[, key]
  #-------------------------------------
  #browser()
  if (fir.st == TRUE) reg <- glm(y ~ ., data = x, family = binomial(link = "logit")) else reg <- glm.nb(y ~ ., data = x)
  y_mat  <- exp(predict(reg, new.data = x_mean)) / (1 + exp(predict(reg, new.data = x_mean)))
  y_is <- y - y_mat + mean(y_mat)
  x_mid <- cbind(y_is, x_non_need)
  for(i in 1:length(x_mid)){
    x.handle <- x_mid[, i]
    if (is.factor(x.handle) == TRUE) x.handle <- as.numeric(x.handle)
    reg2 <-  lm((2 * var.R / mean(x.handle, na.rm = TRUE)) * x.handle ~ R.inc)
    res1[i,1] <- mean(x.handle, na.rm = TRUE)
    res1[i,2] <- sd(x.handle, na.rm = TRUE)
    res1[i,3:4] <- quantile(x.handle, probs = c(0.25, 0.75), na.rm = TRUE)
    res2 <- rbind(res2,cbind(t(summary(reg2)$coefficients[2, c(1, 4)]), t(confint(reg2)[2, ])))
  }
  res_sep <- cbind(res1, res2)
  #------------------
  for(i in 1:length(x_mid)){
    if (res_sep[i,6] <= 0.001) res_sep$sig[i] <- "***"
    else
      if(res_sep[i,6] <= 0.01) res_sep$sig[i] <- "**"
      else
        if (res_sep[i,6] <= 0.05) res_sep$sig[i] <- "*"
        else
          if (res_sep[i,6] <= 0.1) res_sep$sig[i] <- "."
          else res_sep$sig[i] <- " "
  }
  #-------------------
  row.names(res_sep) <- c("y_HI", names(x_non_need))
  names(res_sep) <- c("mean", "sd", "x2.5 %", "x97.5 %", "Estimate", "Pr(>|t|)", "ci2.5 %", "ci97.5 %", "sig")
  #print(res_sep)
  # this part is for res_com
  res_com <- data.frame("mean(sd)" = character("0"), "CI" = character("0"))
  res_com$mean.sd. <- as.character(res_com$mean.sd.)
  res_com$CI <- as.character(res_com$CI)
  for (i in 1:length(x_mid)){
    res_com[i, 1] <- str_c(round(res_sep[i, 1], 4), "(", round(res_sep[i, 2], 4), ")")
    res_com[i, 2] <- str_c(round(res_sep[i, 5], 4), "(", round(res_sep[i, 7], 4), ", ", round(res_sep[i, 8], 4), ")")
  }
  row.names(res_com) <- row.names(res_sep)
  res_sep <- format(res_sep, justify = "left", digits = 4)
  #this part is for regression
  #reg3 <- lm(y_is ~ ., data = x_non_need)
  #this part is for decomposition
  res3 <- fun_ci_sd2(inc = inc, y = y, x = x, fir.st = fir.st, nam_adj = nam_adj)

  decom <- res3$decom[c(key, TRUE), ]
  #decom <- as.data.frame(summary(reg3)$coefficients[-1, ])
  #decom$contri.coef <- decom$Estimate * as.numeric(res_sep$mean[-1]) / mean(y_is, na.rm = TRUE)
  #decom$contri <- decom$contri.coef * as.numeric(res_sep$Estimate[-1])
  #decom <- rbind(decom, c(0, 0, 0, NA, 0, 0))
  #decom$contri[length(x_non_need) + 1] <- as.numeric(res_sep$Estimate[1]) - sum(decom$contri, na.rm = TRUE)
  decom$contri.rate <- as.numeric(decom$contri) * 100 / abs(as.numeric(res_sep$Estimate[1])) # this rate is based on HI
  #row.names(decom)[length(x_non_need) + 1] <- "residual"
  #-------------
  #decom$sig <- car::recode(as.numeric(decom$`P>|z|`), "0:0.001 = '***'; 0.001:0.01 = '**'; 0.01:0.05 = '*'; 0.05:0.1 = '.'; else = ' '")
  decom <- format(decom, justify = "left", digits = 4)
  res <- list(res_sep, res_com, decom)
  names(res) <- c("res_sep", "res_com", "decom")
  return(res)
}
