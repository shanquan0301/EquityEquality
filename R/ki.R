#' @title KaKwani Index (KI) and Decomposition

#' @description Calculate the KaKwani index, and decompose it among factors.

#' @param inc Variable based to rank y.
#' @param y is the dependent variable.
#' @param x is the data frame of independent variables.
#' @param fir.st if whether fit for first step. Default is TRUE.
#' @param nam_adj variables needed to be transfered into dummy variavles. Variables list in nam_adj must be factors.
#' @param need_fa variables used to reflect the health need.


#' @author Shanquan CHEN \email{shanquan0301@gmial.com}

#' @keywords KaKwani Index
#' @examples
##' res <- dat_2011_rural_imp %$%
##'             fun_ki_sd2(inc = need_comp,
##'             y = inp_yes,
##'             x = data.frame(health_professionals.10.000., te,
##'                          edu, live_alone, nation, smoke_year, drink_year, phy_act, liv_con, mental_num, fam_n_inc,
##'                          hous_inc, ncmi, pens, comi,
##'                          age_dec, gender, self_stat, adl_heal, pain_num, chro_num),
##'             nam_adj = c("edu", "self_stat"),
##'             need_fa = c("age_dec", "gender", "self_stat", "adl_heal", "pain_num", "chro_num"))
##'
##'
##'

#' @import dplyr magrittr stringr stats mfx MASS car
#' @importFrom tidyr spread
#' @importFrom nnet class.ind

#' @export fun_ki_sd2
#'
#----------------------------------------------
fun_ki_sd2 <- function(inc, y, x, fir.st = TRUE, nam_adj = nam_adj, need_fa = need_fa){
  #------------------------------
  rank.inc.kna <- rank(inc, na.last = "keep")
  rank.inc.rna <- rank(inc, na.last = NA)
  R.inc <- (rank.inc.kna - 0.5) / length(rank.inc.rna)
  var.R <- var(R.inc, na.rm = TRUE)
  #-------------------------------
  #---------
  for(i in 1: length(x)){
    if(sum(!x[, i] == 0) == 0)
      x[, i] <- 0.001
  }
  #----------
  # sub the non need factor in x_dam
  key <- NA
  for (i in 1:length(names(x))){
    if (sum(str_detect(names(x)[i], need_fa)) > 0)
      key[i] <- FALSE
    else
      key[i] <- TRUE
  }
  x_non_need <- x[, key]
  #---------------
  key <- NA
  for (i in 1:length(nam_adj)){
    if (sum(str_detect(nam_adj[i], need_fa)) > 0)
      key[i] <- FALSE
    else
      key[i] <- TRUE
  }
  #----------------------------------
  #calculate the CI and Gini separately
  res1 <- fun_ci_sd2(inc, y, x_non_need, fir.st = fir.st, nam_adj = nam_adj[key])
  row.names(res1$res_sep)[1] <- "y_CI"
  row.names(res1$res_com)[1] <- "y_CI"
  res2 <- fun_ci_sd3(inc, inc, x_non_need, fir.st = fir.st, nam_adj = nam_adj[key])
  row.names(res2$res_sep)[1] <- "Gini"
  row.names(res2$res_com)[1] <- "Gini"
  #----------------------------------
  #calculate the KI
  reg <-  lm((2 * var.R * (y/mean(y) - inc/mean(inc)))  ~ R.inc)
  #----------------------------------
  #this part is for decomposition
  decom <- cbind(res1$res_sep, res1$res_com[, 1:2])
  decom <- rbind(decom, c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
  decom$Estimate[1] <- reg$coefficients[2]
  decom$`Pr(>|t|)`[1] <- summary(reg)$coefficients[2, 4]
  decom$`ci2.5 %`[1] <- confint(reg)[2, 1]
  decom$`ci97.5 %`[1] <- confint(reg)[2, 2]
  decom$CI[1] <- str_c(round(as.numeric(decom[1, 5]), 4), "(", round(as.numeric(decom[1, 7]), 4), ", ", round(as.numeric(decom[1, 8]), 4), ")")
  decom$sig[1] <- car::recode(as.numeric(decom$`Pr(>|t|)`[1]), "0:0.001 = '***'; 0.001:0.01 = '**'; 0.01:0.05 = '*'; 0.05:0.1 = '.'; else = ' '")
  decom$contri.coef <- c(NA, as.numeric(res1$decom$contri.coef) - as.numeric(res2$decom$contri.coef))
  decom$contri <- c(NA, as.numeric(res1$decom$contri) - as.numeric(res2$decom$contri))
  decom$contri.rate <- decom$contri * 100 / abs(as.numeric(decom$Estimate[1])) # this rate is based on KI
  row.names(decom)[1] <- "y_KI"
  row.names(decom)[length(fun_x_dam2(x_non_need, nam_adj = nam_adj[key])) + 2] <- "residual"
  #-------------
  decom$Estimate <- as.numeric(decom$Estimate)
  decom$`Pr(>|t|)` <- as.numeric(decom$`Pr(>|t|)`)
  decom$`ci2.5 %` <- as.numeric(decom$`ci2.5 %`)
  decom$`ci97.5 %` <- as.numeric(decom$`ci97.5 %`)
  decom <- format(decom, justify = "left", digits = 4)
  res <- list(res1, res2, decom)
  names(res) <- c("res_CI", "res_Gini", "decom")
  return(res)
}
#----------------------------
#only decom the CI
fun_ki_sd2.1 <- function(inc, y, x, fir.st = TRUE, nam_adj = nam_adj, need_fa = need_fa){
  #------------------------------
  rank.inc.kna <- rank(inc, na.last = "keep")
  rank.inc.rna <- rank(inc, na.last = NA)
  R.inc <- (rank.inc.kna - 0.5) / length(rank.inc.rna)
  var.R <- var(R.inc, na.rm = TRUE)
  #---------
  for(i in 1: length(x)){
    if(sum(!x[, i] == 0) == 0)
      x[, i] <- 0.001
  }
  #----------
  #-------------------------------
  # sub the non need factor in x
  key <- NA
  for (i in 1:length(names(x))){
    if (sum(str_detect(names(x)[i], need_fa)) > 0)
      key[i] <- FALSE
    else
      key[i] <- TRUE
  }
  x_non_need <- x[, key]
  #---------------
  key <- NA
  for (i in 1:length(nam_adj)){
    if (sum(str_detect(nam_adj[i], need_fa)) > 0)
      key[i] <- FALSE
    else
      key[i] <- TRUE
  }
  #----------------------------------
  #calculate the CI and Gini separately
  res1 <- fun_ci_sd2(inc, y, x_non_need, fir.st = fir.st, nam_adj = nam_adj[key])
  row.names(res1$res_sep)[1] <- "y_CI"
  row.names(res1$res_com)[1] <- "y_CI"
  res2 <- fun_ci_sd3(inc, inc, x_non_need, fir.st = fir.st, nam_adj = nam_adj[key])
  row.names(res2$res_sep)[1] <- "Gini"
  row.names(res2$res_com)[1] <- "Gini"
  #----------------------------------
  #calculate the KI
  reg <-  lm((2 * var.R * (y/mean(y) - inc/mean(inc)))  ~ R.inc)
  #----------------------------------
  #this part is for decomposition
  decom <- cbind(res1$res_sep, res1$res_com[, 1:2])
  decom <- rbind(decom, c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
  decom$Estimate[1] <- reg$coefficients[2] # this estimate is KI but not the coefficients of regressoiom
  decom$`Pr(>|t|)` <- NA
  decom$`Pr(>|t|)` <- c(NA, res1$decom$`P>|z|`) # let the p inhert from decomposition of CI_need
  #decom$`Pr(>|t|)`[1] <- summary(reg)$coefficients[2, 4]
  decom$`ci2.5 %`[1] <- confint(reg)[2, 1]
  decom$`ci97.5 %`[1] <- confint(reg)[2, 2]
  decom$CI[1] <- str_c(round(as.numeric(decom[1, 5]), 4), "(", round(as.numeric(decom[1, 7]), 4), ", ", round(as.numeric(decom[1, 8]), 4), ")")
  decom$sig <- car::recode(as.numeric(decom$`Pr(>|t|)`), "0:0.001 = '***'; 0.001:0.01 = '**'; 0.01:0.05 = '*'; 0.05:0.1 = '.'; else = ' '")
  decom$contri.coef <- c(NA, as.numeric(res1$decom$contri.coef))
  decom$contri.coef[which(row.names(decom) == "need_comp")] <-  decom$contri.coef[which(row.names(decom) == "need_comp")] - 1
  decom$contri <- c(NA, as.numeric(res1$decom$contri))
  decom$contri[which(row.names(decom) == "need_comp")] <- decom$contri[which(row.names(decom) == "need_comp")] - as.numeric(res2$res_sep$Estimate[1])
  decom$contri.rate <- decom$contri * 100 / abs(as.numeric(decom$Estimate[1])) # this rate is based on KI
  row.names(decom)[1] <- "y_KI"
  row.names(decom)[length(fun_x_dam2(x_non_need, nam_adj = nam_adj[key])) + 2] <- "residual"
  #-------------
  decom$Estimate <- as.numeric(decom$Estimate)
  decom$`Pr(>|t|)` <- as.numeric(decom$`Pr(>|t|)`)
  decom$`ci2.5 %` <- as.numeric(decom$`ci2.5 %`)
  decom$`ci97.5 %` <- as.numeric(decom$`ci97.5 %`)
  decom <- format(decom, justify = "left", digits = 4)
  res <- list(res1, res2, decom)
  names(res) <- c("res_CI", "res_Gini", "decom")
  return(res)
}
