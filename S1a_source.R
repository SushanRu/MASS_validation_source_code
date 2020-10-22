#######################################
############ functions for testing models ##############
#######################################

# a function for significance test, df.dif is the difference between df of model.test and model.ori, defalt = 1
# result is significant sig < 0.1 
sig <- function(model.test, model.ori, df.dif = 1){
  sig <- 1 - pchisq(2*(model.test$logl-model.ori$logl),df.dif)
  result <- c(sig, model.test$logl - model.ori$logl, model.test$logl, model.ori$logl)
  names(result) <- c("sig", "logl.diff", "m.t.logl", "m.ori.logl")
  return(result)
}

###########################################################################
######   function to calculate marker effects for each individual  ########
###########################################################################
# input: H1.geno, H2.geno, marker.name
# output: H1.effect, H2.effect
marker.effect.function <- function(H1.geno, H2.geno, marker.name, Coef.marker, markerEffect){ # maker.name: Ma
  marker1 <- paste(marker.name, "H1", sep = "") # e.g., MaH1
  marker2 <- paste(marker.name, "H2", sep = "") # e.g., MaH2
  H1.effect <- H2.effect <- rep(0, length(H1.geno))
  for (allele in as.character(markerEffect[, marker1])){
    H1.effect[which(H1.geno == allele)] <- Coef.marker[paste(marker1, "_" ,allele, sep = "")]
    H2.effect[which(H2.geno == allele)] <- Coef.marker[paste(marker1, "_", allele, sep = "")]
  }
  return(list("H1.effect" = H1.effect, "H2.effect" = H2.effect))
}

###########################################
######    function to summarize TA    #####
###########################################
# input: Mbest (model), adata
# output: TA.summary
TA.summarize.func <- function(Mbest, adata){
  data_frmt <- adata
  print(popmean <- Mbest$coefficients$fixed['(Intercept)']) #(Intercept) 24.72396
  Year2010 <- Mbest$coefficients$fixed[2]
  # YearPl_2010 YearPl_2011 (Intercept) 
  # 0.000000    4.777417   24.723962
  
  # fixed effect
  YearPl <- data_frmt$YearPl
  Year.effect <- rep(0, length(YearPl)) 
  Year.effect[which(YearPl == "2011")] <- Year2010
  
  print(var.Mbest <- summary(Mbest)$var)
  # gamma component std.error  z.ratio constraint
  # ped(ID)!ped   0.72045523 34.411424  19.77179 1.740431   Positive
  # MaH1!MaH1.var 0.09578577  4.575058   3.83682 1.192409   Positive
  # R!variance    1.00000000 47.763446  10.80735 4.419534   Positive
  print(VM <- var.Mbest$component[2]) # 4.575058
  print(VA <- var.Mbest$component[1] + var.Mbest$component[2]) # 38.98648
  print(VE <- var.Mbest$component[3]) # 47.76345
  print(h2 <- VA/(VA+VE)) # 0.4494123
  print(p <- VM/VA) # 0.1173499
  # variance component
  v.p <- list("var" = var.Mbest, "VM"= VM, "VA" = VA, "VE" = VE, "h2" = h2, "p" = p) 
  
  # marker effects
  Mbest.pvs <- predict(Mbest, classify="MaH1", only = "MaH1")
  markerEffect <- Mbest.pvs$predictions$pvals # prediction of allele effect (for random effect)
  # predict background additive effect
  backg <- predict(Mbest, classify = "ID") # prediction of background additive effect
  ID.predic <- as.data.frame(backg$predictions$pvals)[-c(1:83),] # remove 83 estimates on ancestors
  Coef.ID <- backg$coefficients$random[-c(1:83)][1:321] # coefficient of ID effect on 1001-1709
  Coef.marker <- backg$coefficients$random[-c(1:83)][-c(1:321)]
  fit <- fitted(Mbest) # = breeding values. phenotype = fitted_value + residue, fitted_value = breeding_Value = intercept + ID_coefficient + PgH1 + PgH2
  resid <- resid(Mbest)
  # marker effects
  MaH1.geno <- levels(data_frmt$MaH1)[data_frmt$MaH1]
  MaH2.geno <- levels(data_frmt$MaH2)[data_frmt$MaH2]
  marker.name <- "Ma"
  temp <- marker.effect.function(MaH1.geno, MaH2.geno, marker.name, Coef.marker, markerEffect)
  MaH1.effect <- temp$H1.effect
  MaH2.effect <- temp$H2.effect
  # summary
  phenotype <- data_frmt$TA10
  pop.mean.temp <- rep(popmean, length(phenotype))
  summary <- cbind(phenotype, "popmean"= pop.mean.temp, Year.effect, Coef.ID, MaH1.effect, MaH2.effect, resid, fit, MaH1.geno, MaH2.geno)
  write.csv(summary, file = "TA.summary.csv")
  
  output <- list("model" = TA.best, "popmean" = popmean, "fixed" = Mbest$coefficients$fixed, "v.p" = v.p, "summary" = summary, "markerEffect" = markerEffect)
  return(output)
}

###########################################
######    function to summarize Cn    #####
###########################################
# input: Mbest (model), adata
# output: Cn.summary
Cn.summarize.func <- function(Mbest, adata){
  data_frmt <- adata
  
  print(popmean <- Mbest$coefficients$fixed) #(Intercept) 123.0417
  
  # heritability
  print(var.Mbest <- summary(Mbest)$var)
  print(VM <- var.Mbest$component[2]) # 161.006
  print(VA <- var.Mbest$component[1] + var.Mbest$component[2]) #988.1767
  print(VE <- var.Mbest$component[3]) #3303.172
  print(h2 <- VA/(VA+VE)) # 0.2303
  print(p <- VM/VA) # 0.1629
  # variance component
  v.p <- list("var" = var.Mbest, "VM"= VM, "VA" = VA, "VE" = VE, "h2" = h2, "p" = p) 
  
  # marker effects
  Mbest.pvs <- predict(Mbest, classify="PgH1", only = "PgH1")
  markerEffect <- Mbest.pvs$predictions$pvals # prediction of allele effect (for random effect)
  
  # predict background additive effect
  backg <- predict(Mbest, classify = "ID") # prediction of background additive effect
  ID.predic <- as.data.frame(backg$predictions$pvals)[-c(1:83),]
  Coef.ID <- backg$coefficients$random[-c(1:83)][1:321] # coefficient of ID effect on 1001-1709
  Coef.marker <- backg$coefficients$random[-c(1:83)][-c(1:321)]
  fit <- fitted(Mbest) # = breeding values, phenotype = fitted_value + residue, fitted_value = breeding_Value = intercept + ID_coefficient + PgH1 + PgH2
  resid <- resid(Mbest)
  # marker effects
  PgH1.geno <- levels(data_frmt$PgH1)[data_frmt$PgH1]
  PgH2.geno <- levels(data_frmt$PgH2)[data_frmt$PgH2]
  marker.name <- "Pg"
  temp <- marker.effect.function(PgH1.geno, PgH2.geno, marker.name, Coef.marker, markerEffect)
  PgH1.effect <- temp$H1.effect
  PgH2.effect <- temp$H2.effect
  
  # summary
  phenotype <- data_frmt$ES_Cn
  pop.mean.temp <- rep(popmean, length(phenotype))
  summary <- cbind(phenotype, "popmean" = pop.mean.temp, Coef.ID, PgH1.effect, PgH2.effect, resid, fit, PgH1.geno, PgH2.geno)
  write.csv(summary, file = "Cn.summary.csv")
  
  output <- list("model" = Cn.best, "popmean" = popmean, "fixed" = Mbest$coefficients$fixed, "v.p" = v.p, "summary" = summary, "markerEffect" = markerEffect)
  return(output)
}

###########################################
######    function to summarize M2    #####
###########################################
# input: Mbest (model), adata
# output: TA.summary
M2.summarize.func <- function(Mbest, adata){
  data_frmt <- adata
  print(popmean <- Mbest$coefficients$fixed['(Intercept)']) #(Intercept) 17.36956
  
  print(var.Mbest <- summary(Mbest)$var)
  # gamma component std.error   z.ratio constraint
  # ped(ID)!ped   0.83853582 7.4823467 3.7622957 1.9887716   Positive
  # MaH1!MaH1.var 0.03408447 0.3041394 0.3603502 0.8440107   Positive
  # PgH1!PgH1.var 0.14093336 1.2575638 1.3913139 0.9038678   Positive
  # R!variance    1.00000000 8.9231092 2.0270965 4.4019163   Positive
  
  print(VM.Ma <- var.Mbest$component[2]) # 0.3041394
  print(VM.Pg <- var.Mbest$component[3]) # 1.257564
  print(VA <- sum(var.Mbest$component[1:3])) # 9.04405
  print(VE <- var.Mbest$component[4]) # 8.923109
  print(h2 <- VA/(VA+VE)) # 0.5033656
  print(p <- (VM.Ma + VM.Pg)/VA) # 0.1726774
  
  # variance component
  v.p <- list("var" = var.Mbest, "VM.Ma"= VM.Ma, "VM.Pg" = VM.Pg, "VA" = VA, "VE" = VE, "h2" = h2, "p" = p) 
  
  # marker effects
  Mbest.pvs.Ma <- predict(Mbest, classify="MaH1", only = "MaH1")
  markerEffect.Ma <- Mbest.pvs.Ma$predictions$pvals # prediction of allele effect (for random effect)
  
  Mbest.pvs.Pg <- predict(Mbest, classify="PgH1", only = "PgH1")
  markerEffect.Pg <- Mbest.pvs.Pg$predictions$pvals # prediction of allele effect (for random effect)
  
  # predict background additive effect
  backg <- predict(Mbest, classify = "ID") # prediction of background additive effect
  ID.predic <- as.data.frame(backg$predictions$pvals)[-c(1:83),] # remove 83 estimates on ancestors
  Coef.ID <- backg$coefficients$random[-c(1:83)][1:321] # coefficient of ID effect on 1001-1709
  Coef.marker <- backg$coefficients$random[-c(1:83)][-c(1:321)]
  fit <- fitted(Mbest) # = breeding values. phenotype = fitted_value + residue, fitted_value = breeding_Value = intercept + ID_coefficient + PgH1 + PgH2
  resid <- resid(Mbest)
  # marker effects
  MaH1.geno <- levels(data_frmt$MaH1)[data_frmt$MaH1]
  MaH2.geno <- levels(data_frmt$MaH2)[data_frmt$MaH2]
  marker.name <- "Ma"
  temp <- marker.effect.function(MaH1.geno, MaH2.geno, marker.name, Coef.marker, markerEffect.Ma)
  MaH1.effect <- temp$H1.effect
  MaH2.effect <- temp$H2.effect
   
  PgH1.geno <- levels(data_frmt$PgH1)[data_frmt$PgH1]
  PgH2.geno <- levels(data_frmt$PgH2)[data_frmt$PgH2]
  marker.name <- "Pg"
  temp <- marker.effect.function(PgH1.geno, PgH2.geno, marker.name, Coef.marker, markerEffect.Pg)
  PgH1.effect <- temp$H1.effect
  PgH2.effect <- temp$H2.effect
  
  # summary
  phenotype <- data_frmt$ES_M2_lb
  pop.mean.temp <- rep(popmean, length(phenotype))
  summary <- cbind(phenotype, "popmean"= pop.mean.temp, Coef.ID, MaH1.effect, MaH2.effect, PgH1.effect, PgH2.effect, resid, fit, MaH1.geno, MaH2.geno, PgH1.geno, PgH2.geno)
  write.csv(summary, file = "Pg.summary.csv")
  
  output <- list("model" = TA.best, "popmean" = popmean, "fixed" = Mbest$coefficients$fixed, "v.p" = v.p, "summary" = summary, "markerEffect.Ma" = markerEffect.Ma,
                 "markerEffect.Pg" = markerEffect.Pg)
  return(output)
}


##########################################################
#######    function to write results into file     #######
##########################################################
# input: TA.summary, M2.summary, Cn.summary, filename
# output: none

write.func <- function(result, filename){
  TA.summary <- result$TA.summary
  M2.summary <- result$M2.summary
  Cn.summary <- result$Cn.summary
  
  cat("Model, TA100~1 + YearPl random = ~ ped(ID) + MaH1 + and(MaH2)\n", file = filename)
  cat("popmean,", TA.summary$popmean, "\n",file = filename, append = T)
  cat("\n", file = filename, append = T)
  cat("fixed effect", file = filename, append = T)
  write.table(TA.summary$fixed, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  cat("marker effect", file = filename, append = T)
  write.table(TA.summary$markerEffect, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  write.table(TA.summary$v.p$var, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  write.table(TA.summary$v.p[-1], file = filename, col.names = NA, sep = ",",append = T)
  
  cat("\n", file = filename, append = T)
  cat("\n", file = filename, append = T)
  cat("Model, ES_M2_lb ~ 1 random = ~ ped(ID) + MaH1 + and(MaH2) + PgH1 + and(PgH2)\n", file = filename, append = T)
  cat("popmean,", M2.summary$popmean, "\n",file = filename, append = T)
  cat("\n", file = filename, append = T)
  cat("fixed efefct", file = filename, append = T)
  write.table(M2.summary$fixed, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  cat("marker effect\n", file = filename, append = T)
  write.table(M2.summary$markerEffect.Ma, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  write.table(M2.summary$markerEffect.Pg, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  write.table(M2.summary$v.p$var, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  write.table(M2.summary$v.p[-1], file = filename, col.names = NA, sep = ",",append = T)
  
  cat("\n", file = filename, append = T)
  cat("\n", file = filename, append = T)
  cat("Model, ES_Cn ~ 1 random = ~ ped(ID) + PgH1+and(PgH2)\n", file = filename, append = T)
  cat("popmean,", Cn.summary$popmean, "\n",file = filename, append = T)
  cat("\n", file = filename, append = T)
  cat("fixed effect", file = filename, append = T)
  write.table(Cn.summary$fixed, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  cat("marker effect", file = filename, append = T)
  write.table(Cn.summary$markerEffect, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  write.table(Cn.summary$v.p$var, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  write.table(Cn.summary$v.p[-1], file = filename, col.names = NA, sep = ",",append = T)
}