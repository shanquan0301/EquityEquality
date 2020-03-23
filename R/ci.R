#' @title Concentration Index (CI) and Decomposition

#' @description Calculate the concentration index, and decompose it among factors.

#' @param inc Variable based to rank y.
#' @param y is the dependent variable.
#' @param x is the data frame of independent variables.
#' @param fir.st if whether fit for first step. Default is TRUE.
#' @param nam_adj variables needed to be transfered into dummy variavles. Variables list in nam_adj must be factors.


#' @author Shanquan CHEN \email{shanquan0301@gmial.com}

#' @keywords Concentration Index
#' @examples
##' res <- dat_2011_rural_imp %$%
##'            fun_ci_sd2(inc = hous_inc,
##'            y = outp_yes,
##'            x = data.frame(health_professionals.10.000., te,
##'                           edu, live_alone, nation, smoke_year, drink_year, phy_act, liv_con, mental_num, fam_n_inc,
##'                           hous_inc, ncmi, pens, comi,
##'                           age_dec, gender, self_stat, adl_heal, pain_num, chro_num),
##'            nam_adj = c("edu", "self_stat"))

##' res <- dat_2011_rural_imp %>% filter(outp_tim > 0) %$%
##'          fun_ci_sd2(inc = hous_inc,
##'          y = outp_tim,
##'          x = data.frame(health_professionals.10.000., te,
##'                         edu, live_alone, nation, smoke_year, drink_year, phy_act, liv_con, mental_num, fam_n_inc,
##'                         hous_inc, ncmi, pens, comi,
##'                         age_dec, gender, self_stat, adl_heal, pain_num, chro_num),
##'             fir.st = FALSE,
##'             nam_adj = c("edu", "self_stat"))



#' @import dplyr magrittr stringr stats mfx MASS car
#' @importFrom tidyr spread
#' @importFrom nnet class.ind

#' @export fun_ci_sd2
#' @export fun_ci_sd3
#'
#----------------------------------------------
#fun_ci_sd2 is for logit and NB
fun_ci_sd2 <- function(inc, y, x, fir.st = TRUE, nam_adj = nam_adj){
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
  x_mid <- cbind(y, x)
  x <- fun_x_dam2(data = x, nam_adj = nam_adj)
  x_mid <- fun_x_dam2(data = x_mid, nam_adj = nam_adj)
  # this part is for res_sep-----------------------------------------------------------
  for(i in 1:length(x_mid)){
    x.handle <- x_mid[, i]
    if (is.factor(x.handle) == TRUE) x.handle <- as.numeric(x.handle)
    reg <-  lm((2 * var.R / mean(x.handle, na.rm = TRUE)) * x.handle ~ R.inc)
    res1[i,1] <- mean(x.handle, na.rm = TRUE)
    res1[i,2] <- sd(x.handle, na.rm = TRUE)
    res1[i,3:4] <- quantile(x.handle, probs = c(0.25, 0.75), na.rm = TRUE)
    res2 <- rbind(res2,cbind(t(summary(reg)$coefficients[2, c(1, 4)]), t(confint(reg)[2, ])))
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
  row.names(res_sep) <- names(x_mid)
  #browser()
  names(res_sep) <- c("mean", "sd", "x2.5 %", "x97.5 %", "Estimate", "Pr(>|t|)", "ci2.5 %", "ci97.5 %", "sig")
  # this part is for res_com and res_table--------------------------------------------------------
  res_com <- data.frame("mean(sd)" = character("0"), "CI" = character("0"))
  res_com$mean.sd. <- as.character(res_com$mean.sd.)
  res_com$CI <- as.character(res_com$CI)
  res_table <- list()
  for (i in 1:length(x_mid)){
    res_com[i, 1] <- str_c(round(res_sep[i, 1], 4), "(", round(res_sep[i, 2], 4), ")")
    res_com[i, 2] <- str_c(round(res_sep[i, 5], 4), "(", round(res_sep[i, 7], 4), ", ", round(res_sep[i, 8], 4), ")")
    res_table[[i]] <- table(x_mid[, i])
  }
  names(res_table) <- names(x_mid)
  res_sep <- format(res_sep, justify = "left", digits = 4)
  #this part is for regression
  #transfer factor x into numeric
  for (i in 1:length(x)){
    if (is.factor(x[, i])) x[, i] <- as.numeric(x[, i])
  }
  #formula <- as.formula(str_c(c("y", str_c(names(x), collapse = " + ")), collapse = " ~ "))
  set.seed(1000)
  #browser()
  if (fir.st == TRUE) reg <- logitmfx(y ~ ., data = x, robust= TRUE) else reg <- negbinmfx(y ~ ., data = x, robust = TRUE)
  if (fir.st == TRUE) reg2 <- glm(y ~ ., data = x, family = binomial(link = "logit")) else reg2 <- glm.nb(y ~ ., data = x)
  #this part is for decomposition
  #browser()
  decom <- as.data.frame(reg$mfxest)
  decom$name <- row.names(decom)
  mdat <- data.frame(name = names(reg2$coefficients)[-1], rank = 1: length(names(reg2$coefficients)[-1]))
  decom <- merge(mdat, decom, by = "name", all.x = TRUE)
  decom <- decom[order(decom$rank),]
  row.names(decom) <- decom$name
  decom <- decom[, -c(1:2)]
  #decom$`P>|z|`[which(is.na(decom$`P>|z|`))] <- 0
  #browser()
  decom$contri.coef <- decom$`dF/dx` * as.numeric(res_sep$mean[-1]) / mean(y, na.rm = TRUE)
  decom$contri <- decom$contri.coef * as.numeric(res_sep$Estimate[-1])
  decom <- rbind(decom, c(0, 0, 0, NA, 0, 0))
  reg3 <-  lm((2 * var.R / mean(y, na.rm = TRUE)) * y ~ R.inc)
  decom$contri[length(x) + 1] <- (reg3$coefficients[2] - sum(decom$contri, na.rm = TRUE))
  decom$contri.rate <- decom$contri * 100 / abs(reg3$coefficients[2]) # this rate is based on CI
  row.names(decom)[length(x) + 1] <- "residual"
  #-------------
  decom$sig <- car::recode(decom$`P>|z|`, "0:0.001 = '***'; 0.001:0.01 = '**'; 0.01:0.05 = '*'; 0.05:0.1 = '.'; else = ' '")

  #------------
  res_reg2 <- as.data.frame(summary(reg2)$coefficients)
  res_reg2$name <- row.names(res_reg2)
  mdat <- data.frame(name = names(reg2$coefficients)[-1], rank = 1: length(names(reg2$coefficients)[-1]))
  res_reg2 <- merge(mdat, res_reg2, by = "name", all.x = TRUE)
  res_reg2 <- res_reg2[order(res_reg2$rank),]
  p.value <- res_reg2$`Pr(>|z|)`
  #p.value[which(is.na(p.value))] <- 0
  p.value <- c(p.value, NA)
  decom$sig_ori <- car::recode(p.value, "0:0.001 = '***'; 0.001:0.01 = '**'; 0.01:0.05 = '*'; 0.05:0.1 = '.'; else = ' '")

  #-----------------
  decom <- format(decom, justify = "left", digits = 4)

  #-------------------
  row.names(res_com) <- names(x_mid)
  names(res_com) <- c("mean(sd)", "CI")
  res <- list(res_sep, res_com, res_table, reg2, reg, decom)
  names(res) <- c("res_sep", "res_com", "res_table", "reg_ori", "reg_mfx", "decom")
  return(res)
  #print("res_sep")
  #print(res_sep, right = FALSE)
  #if (table == TRUE) print("res_com")
  #print(res_com)
}

#fun_ci_sd3 is for linear regression
fun_ci_sd3 <- function(inc, y, x, fir.st = TRUE, nam_adj = nam_adj){
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
  x_mid <- cbind(y, x)
  x <- fun_x_dam2(data = x, nam_adj = nam_adj)
  x_mid <- fun_x_dam2(data = x_mid, nam_adj = nam_adj)
  # this part is for res_sep-----------------------------------------------------------
  #if (is.data.frame(x) == TRUE){
  for(i in 1:length(x_mid)){
    x.handle <- x_mid[, i]
    if (is.factor(x.handle) == TRUE) x.handle <- as.numeric(x.handle)
    reg <-  lm((2 * var.R / mean(x.handle, na.rm = TRUE)) * x.handle ~ R.inc)
    res1[i,1] <- mean(x.handle, na.rm = TRUE)
    res1[i,2] <- sd(x.handle, na.rm = TRUE)
    res1[i,3:4] <- quantile(x.handle, probs = c(0.25, 0.75), na.rm = TRUE)
    res2 <- rbind(res2,cbind(t(summary(reg)$coefficients[2, c(1, 4)]), t(confint(reg)[2, ])))
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
  row.names(res_sep) <- c("y", names(x))
  #browser()
  names(res_sep) <- c("mean", "sd", "x2.5 %", "x97.5 %", "Estimate", "Pr(>|t|)", "ci2.5 %", "ci97.5 %", "sig")
  # this part is for res_com and res_table--------------------------------------------------------
  res_com <- data.frame("mean(sd)" = character("0"), "CI" = character("0"))
  res_com$mean.sd. <- as.character(res_com$mean.sd.)
  res_com$CI <- as.character(res_com$CI)
  res_table <- list()
  for (i in 1:length(x_mid)){
    res_com[i, 1] <- str_c(round(res_sep[i, 1], 4), "(", round(res_sep[i, 2], 4), ")")
    res_com[i, 2] <- str_c(round(res_sep[i, 5], 4), "(", round(res_sep[i, 7], 4), ", ", round(res_sep[i, 8], 4), ")")
    res_table[[i]] <- table(x_mid[, i])
  }
  names(res_table) <- names(x_mid)
  row.names(res_com) <- row.names(res_sep)
  res_sep <- format(res_sep, justify = "left", digits = 4)
  #this part is for regression
  reg3 <- lm(y ~ ., data = x)
  #this part is for decomposition
  #browser()

  decom <- as.data.frame(summary(reg3)$coefficients[-1, ])
  decom$name <- str_remove_all(row.names(decom), "`")
  mdat <- data.frame(name = names(x), rank = 1: length(names(x)))
  decom <- merge(mdat, decom, by = "name", all.x = TRUE)
  decom <- decom[order(decom$rank),]
  row.names(decom) <- decom$name
  decom <- decom[, -c(1:2)]

  decom$contri.coef <- decom$Estimate * as.numeric(res_sep$mean[-1]) / mean(y_is, na.rm = TRUE)
  decom$contri <- decom$contri.coef * as.numeric(res_sep$Estimate[-1])
  decom <- rbind(decom, c(0, 0, 0, NA, 0, 0))
  decom$contri[length(x) + 1] <- as.numeric(res_sep$Estimate[1]) - sum(decom$contri, na.rm = TRUE)
  decom$contri.rate <- decom$contri * 100 / abs(as.numeric(res_sep$Estimate[1])) # this rate is based on CI
  row.names(decom)[length(x) + 1] <- "residual"

  #-------------
  decom$sig <- car::recode(decom$`P>|z|`, "0:0.001 = '***'; 0.001:0.01 = '**'; 0.01:0.05 = '*'; 0.05:0.1 = '.'; else = ' '")
  #-----------------
  decom <- format(decom, justify = "left", digits = 4)
  #-------------------
  res <- list(res_sep, res_com, decom)
  names(res) <- c("res_sep", "res_com", "decom")
  return(res)
}


#---------------------------------------
# inner use functions
# fun_x_dam is to transfer factor into dummy variable
fun_x_dam <- function(data, nam_adj){
  for (i in nam_adj){
    nam <- levels(data[, i])
    data[, nam[-1]] <- 0
    for (j in 1: dim(data)[1]){
      if (data[j, i] != nam[1] & !is.na(data[j, i]))
        data[j, as.character(data[j, i])] <- 1
    }
    loc1 <- which(names(data) == i) - 1
    loc2 <- (length(data) - length(nam) + 2)
    loc3 <- which(names(data) == i) + 1
    loc4 <- loc2 - 1
    nam <- str_c("`", nam)
    nam <- str_c(nam, "`")
    if (loc1 > 0){
      names(data)[loc2:length(data)] <- str_c(i, nam[-1], sep = ":")
      data <- data[, c(1:loc1,  loc2:length(data), loc3:loc4)]
    }
    else{
      names(data)[loc2:length(data)] <- str_c(i, nam[-1], sep = ":")
      data <- data[, c(loc2:length(data), loc3:loc4)]
    }
  }
  del <- NULL
  for (i in 1: length(data)){
    if (sum(as.numeric(data[, i]), na.rm = TRUE) == 0)
      del <- c(del, i)
  }
  if (length(del) > 0) data <- data[, -del]
  return(data)
}

fun_x_dam2 <- function(data, nam_adj){
  for (i in nam_adj){
    loc1 <- which(names(data) == i )
    if(loc1 == 1){
      mdat1 <- data[, names(data)[1]]
      mdat1 <- as.data.frame(class.ind(mdat1))[, -1]
      names(mdat1) <- str_c(i, ": ", colnames(mdat1))
      mdat3 <- data[, names(data)[(loc1 + 1) : (length(names(data)))]]
      data <- cbind(mdat1, mdat3)
    } else {
      mdat1 <- data[, names(data)[1:(loc1-1)]]
      mdat2 <- data[, names(data)[loc1]]
      mdat3 <- data %>% dplyr::select(names(data)[(loc1 + 1) : (length(names(data)))])
      mdat2 <- as.data.frame(class.ind(mdat2))[, -1]
      names(mdat2) <- str_c(i, ": ", colnames(mdat2))
      data <- cbind(mdat1, mdat2, mdat3)
    }
  }
  return(data)
}

