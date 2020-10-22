#########
# functions to analyze marker allele effects and breeding values based on selected models

Mbest.calc <- function(adata, Mbest, random = T){
  popmean <- Mbest$coefficients$fixed["(Intercept)"]
  
  #variance compoent
  if (random == T){
    var.Mbest <- summary(Mbest)$var
    # gamma component std.error   z.ratio constraint
    # ped(ID)!ped   0.83853582 7.4823467 3.7622957 1.9887716   Positive
    # MaH1!MaH1.var 0.03408447 0.3041394 0.3603502 0.8440107   Positive
    # PgH1!PgH1.var 0.14093336 1.2575638 1.3913139 0.9038678   Positive
    # R!variance    1.00000000 8.9231092 2.0270965 4.4019163   Positive
    
    VM.Ma <- var.Mbest$component[2] # 0.3041394
    VM.Pg <- var.Mbest$component[3] # 1.257564
    VA <- sum(var.Mbest$component[1:3]) # 9.04405
    VE <- var.Mbest$component[4] # 8.923109
    h2 <- VA/(VA+VE) # 0.5033656
    p <- (VM.Ma + VM.Pg)/VA # 0.1726774
    v.p = list("var" = var.Mbest, "VM.Ma"= VM.Ma, "VM.Pg" = VM.Pg, "VA" = VA, "VE" = VE, "h2" = h2, "p" = p) 
  } else{
    v.p = NA
  }
  
  # predicted marker performance
  predict.Pg <- predict(Mbest, classify="PgH1", only = "PgH1")
  markerEffect.Pg <- predict.Pg$predictions$pvals # prediction of allele effect (for random effect)
  
  predict.Ma <- predict(Mbest, classify="MaH1", only = "MaH1")
  markerEffect.Ma <- predict.Ma$predictions$pvals
  
  # predict pedigree effect
  backg <- predict(Mbest, classify = "ID") # prediction of background additive effect
  ID.predic <- as.data.frame(backg$predictions$pvals)[-c(1:83),]
  
  Coef.ID <- backg$coefficients$random[-c(1:83)][1:321] # coefficient of ID effect on 1001-1709
  if(random == F){
    Coef.marker <- Mbest$coefficients$fixed[c(1:10)]
    filename <- "M.fix.summary.csv"
  }else{
    Coef.marker <- backg$coefficients$random[-c(1:83)][-c(1:321)]
    filename <- "M.random.summary.csv"
  }
  names(Coef.marker) <- gsub("PgH1_", "", names(Coef.marker))
  names(Coef.marker) <- gsub("MaH1_", "", names(Coef.marker))
  
  PgH1 <- levels(adata$PgH1)[adata$PgH1]
  PgH2 <- levels(adata$PgH2)[adata$PgH2]
  Coef.PgH2 <- Coef.PgH1 <- rep(0, length(PgH1))
  Coef.PgH1 <- Coef.marker[PgH1]
  Coef.PgH2 <- Coef.marker[PgH2]
  
  MaH1 <- levels(adata$MaH1)[adata$MaH1]
  MaH2 <- levels(adata$MaH2)[adata$MaH2]
  Coef.MaH2 <- Coef.MaH1 <- rep(0, length(MaH1))
  Coef.MaH1 <- Coef.marker[MaH1]
  Coef.MaH2 <- Coef.marker[MaH2]
  
  # fitted values
  fit <- fitted(Mbest) # phenotype = fitted_value + residue, fitted_value = breeding_Value = intercept + ID_coefficient + PgH1 + PgH2 + MaH1 + MaH2
  resid <- resid(Mbest)
  
  phenotype <- adata$ES_M2_lb
  summary <- data.frame(phenotype, rep(popmean, length(phenotype)), Coef.ID, Coef.PgH1, Coef.PgH2, Coef.MaH1, Coef.MaH2, resid, fit, PgH1, PgH2, MaH1, MaH2)
  Mbest.result <- list("model" = Mbest, "popmean" = popmean, "fixed" = Mbest$coefficients$fixed, "v.p" = v.p, "marker.effect.Pg" = markerEffect.Pg, "marker.effect.Ma" = markerEffect.Ma, "summary" = summary)
  
  write.csv(summary, file = filename)
  write.table(markerEffect.Ma, file = filename, sep = ",", col.names = NA, append = T)
  write.table (markerEffect.Pg, file = filename, sep = ",", col.names = NA, append = T)
  return(Mbest.result)
}


#####################
###  Cnbest.calc  ###
#####################
Cnbest.calc <- function(adata, Cnbest, random = T){
  Mbest = Cnbest
  popmean <- Mbest$coefficients$fixed["(Intercept)"]
  
  if(random == T){
    var.Mbest <- summary(Mbest)$var
    VM <- var.Mbest$component[2] # 161.006
    VA <- var.Mbest$component[1] + var.Mbest$component[2] #988.1767
    VE <- var.Mbest$component[3] #3303.172
    h2 <- VA/(VA+VE) # 0.2303
    p <- VM/VA # 0.1629
    # variance component
    v.p <- list("var" = var.Mbest, "VM"= VM, "VA" = VA, "VE" = VE, "h2" = h2, "p" = p) 
  }else{
    v.p = NA
  }
  
  
  # predicted marker performance
  predict.Pg <- predict(Mbest, classify="PgH1", only = "PgH1")
  markerEffect.Pg <- predict.Pg$predictions$pvals # prediction of allele effect (for random effect)
  
  # predict background additive effect
  backg <- predict(Mbest, classify = "ID") # prediction of background additive effect
  ID.predic <- as.data.frame(backg$predictions$pvals)[-c(1:83),]
  # write.csv(ID.predic, "Ped_predict.csv")
  
  Coef.ID <- backg$coefficients$random[-c(1:83)][1:321] # coefficient of ID effect on 1001-1709
  if(random == F){
    Coef.marker <- Mbest$coefficients$fixed[c(1:3)]
    filename = "Cn.fix.summary.csv"
  }else{
    Coef.marker <- backg$coefficients$random[-c(1:83)][-c(1:321)]
    filename = "Cn.random.summary.csv"
  }
  names(Coef.marker) <- gsub("PgH1_", "", names(Coef.marker))
  PgH1 <- levels(adata$PgH1)[adata$PgH1]
  PgH2 <- levels(adata$PgH2)[adata$PgH2]
  
  Coef.PgH2 <- Coef.PgH1 <- rep(0, length(PgH1))
  Coef.PgH1 <- Coef.marker[PgH1]
  Coef.PgH2 <- Coef.marker[PgH2]
  
  # save fitted values
  fit <- fitted(Mbest) # phenotype = fitted_value + residue, fitted_value = breeding_Value = intercept + ID_coefficient + PgH1 + PgH2
  resid <- resid(Mbest)
  phenotype <- adata$ES_Cn
  summary <- data.frame(phenotype, rep(popmean, length(phenotype)), Coef.ID, Coef.PgH1, Coef.PgH2, resid, fit, PgH1, PgH2)
  
  Cn.result <- list("model" = Mbest, "popmean" = popmean, "fixed" = Mbest$coefficients$fixed,  "v.p" = v.p ,"marker.effect.Pg" = markerEffect.Pg, "summary" = summary)
  
  write.csv(summary, file = filename)
  write.table (markerEffect.Pg, file = filename, sep = ",", col.names = NA, append = T)
  return(Cn.result)
  
}


TAbest.calc <- function(adata, TAbest, random = T) # random = T if marker effects are random, F if fixed
{
  Mbest <- TAbest
  popmean <- Mbest$coefficients$fixed["(Intercept)"] # 2.472396
  Coef.fix <- Mbest$coefficients$fixed[c("YearPl_2010", "YearPl_2011")] # Coefficient of yearPl
  names(Coef.fix) <- gsub("YearPl_", "", names(Coef.fix))

  if(random == T){
    # variance component
    var.Mbest <- summary(Mbest)$var
    # gamma  component std.error  z.ratio constraint
    # ped(ID)!ped   0.72045523 0.34411424 0.1977179 1.740431   Positive
    # MaH1!MaH1.var 0.09578577 0.04575058 0.0383682 1.192409   Positive
    # R!variance    1.00000000 0.47763446 0.1080735 4.419534   Positive
    VM <- var.Mbest$component[2] # 0.04575058
    VA <- var.Mbest$component[1] + var.Mbest$component[2] # 0.3898648
    VE <- var.Mbest$component[3] # 0.4776345
    h2 <- VA/(VA+VE) # 0.4494123
    p <- VM/VA # 0.1173499
    # variance component
    v.p <- list("var" = var.Mbest, "VM"= VM, "VA" = VA, "VE" = VE, "h2" = h2, "p" = p) 
  }else{
    v.p = NA
  }
  
  
  # predicted marker performance
  predict.Ma <- predict(Mbest, classify="MaH1", only = "MaH1")
  markerEffect.Ma <- predict.Ma$predictions$pvals # prediction of allele effect (for random effect)
  
  # predict background additive effect
  backg <- predict(Mbest, classify = "ID") # prediction of background additive effect
  ID.predic <- as.data.frame(backg$predictions$pvals)[-c(1:83),]
  # write.csv(ID.predic, "Ped_predict.csv")
  
  
  if (random == F){
    Coef.marker <- Mbest$coefficients$fixed[c(1:7)]
    filename <- "TA.fix.summary.csv"
  }else{
    Coef.marker <- backg$coefficients$random[-c(1:83)][-c(1:321)]
    filename <- "TA.random.summary.csv"
  }
  names(Coef.marker) <- gsub("MaH1_", "", names(Coef.marker))
  # calculate marker & year coefficient for each individual
  MaH1 <- levels(adata$MaH1)[adata$MaH1]
  MaH2 <- levels(adata$MaH2)[adata$MaH2]
  yearPl <- levels(adata$YearPl)[adata$YearPl]
  Coef.yearPl <- Coef.MaH2 <- Coef.MaH1 <- rep(0, length(MaH1))
  Coef.MaH1 <- Coef.marker[MaH1]
  Coef.MaH2 <- Coef.marker[MaH2]
  Coef.yearPl <- Coef.fix[yearPl]
  
  Coef.ID <- backg$coefficients$random[-c(1:83)][1:321] # coefficient of ID effect on 1001-1709
  # save fitted values
  fit <- fitted(Mbest) # phenotype = fitted_value + residue, fitted_value = breeding_Value = intercept + ID_coefficient + PgH1 + PgH2
  resid <- resid(Mbest)
  
  phenotype <- adata$TA100
  summary <- data.frame(phenotype, rep(popmean, length(phenotype)), 
                        Coef.ID = Coef.ID, Coef.yearPl, Coef.MaH1, Coef.MaH2, resid, fit, yearPl, MaH1, MaH2)
  
  TA.result <- list("model" = Mbest, "popmean" = popmean, "fixed" = Mbest$coefficients$fixed, "v.p" = v.p, "marker.effect.Ma" = markerEffect.Ma, "summary" = summary)
  
  write.csv(summary, file = filename)
  write.table (markerEffect.Ma, file = filename, sep = ",", col.names = NA, append = T)
  return(TA.result)
}


write.func <- function(predictions, filename){
  TA.summary <- predictions$TA.r
  M2.summary <- predictions$M.r
  Cn.summary <- predictions$Cn.r
  
  cat("Model, TA10~1 + YearPl random = ~ ped(ID) + MaH1 + and(MaH2)\n", file = filename)
  cat("popmean,", TA.summary$popmean, "\n",file = filename, append = T)
  cat("\n", file = filename, append = T)
  cat("fixed effect", file = filename, append = T)
  write.table(TA.summary$fixed, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  cat("marker effect", file = filename, append = T)
  write.table(TA.summary$marker.effect.Ma, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  write.table(TA.summary$v.p$var, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  write.table(TA.summary$v.p[-1], file = filename, col.names = NA, sep = ",",append = T)
  
  
  cat("\n", file = filename, append = T)
  cat("\n", file = filename, append = T)
  cat("Model, ES_Cn ~ 1 random = ~ ped(ID) + PgH1+and(PgH2)\n", file = filename, append = T)
  cat("popmean,", Cn.summary$popmean, "\n",file = filename, append = T)
  cat("\n", file = filename, append = T)
  cat("fixed effect", file = filename, append = T)
  write.table(Cn.summary$fixed, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  cat("marker effect", file = filename, append = T)
  write.table(Cn.summary$marker.effect, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  write.table(Cn.summary$v.p$var, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  write.table(Cn.summary$v.p[-1], file = filename, col.names = NA, sep = ",",append = T)
  
  cat("\n", file = filename, append = T)
  cat("\n", file = filename, append = T)
  cat("Model, ES_M2_lb ~ 1 random = ~ ped(ID) + MaH1 + and(MaH2) + PgH1 + and(PgH2)\n", file = filename, append = T)
  cat("popmean,", M2.summary$popmean, "\n",file = filename, append = T)
  cat("\n", file = filename, append = T)
  cat("fixed efefct", file = filename, append = T)
  write.table(M2.summary$fixed, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  cat("marker effect\n", file = filename, append = T)
  write.table(M2.summary$marker.effect.Ma, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  write.table(M2.summary$marker.effect.Pg, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  write.table(M2.summary$v.p$var, file = filename, col.names = NA, sep = ",",append = T)
  cat("\n", file = filename, append = T)
  write.table(M2.summary$v.p[-1], file = filename, col.names = NA, sep = ",",append = T)
}

