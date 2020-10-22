library(lattice)
library(asreml)
source(file = "S1b_source.R")
######################################################################
# S2. estimate marker and genotypeic effect based on selected models #
######################################################################
# load data
load(file = "data.frmt.rda") # input data
adata <- input$data_frmt
load(file = "Model.select.rda", verbose = T)
TAbest <- Model.select$TAbest
Cnbest <- Model.select$Cnbest
Mbest <- Model.select$Mbest
TA.fix <- Model.select$TA.fix
Cn.fix <- Model.select$Cn.fix
M.fix <- Model.select$M.fix
  
# # best models

# 
# Cnbest <- Cn08 <- asreml(ES_Cn ~ 1,
#                random = ~ ped(ID) + PgH1+and(PgH2),
#                ginverse=list(ID=Ainv),
#                data=adata)
# 
# Mbest <- asreml(ES_M2_lb ~ 1,
#               random = ~ ped(ID) + MaH1 + and(MaH2) + PgH1 + and(PgH2),
#               ginverse=list(ID=Ainv),
#               data=adata)
# 
# # fixed effect models
# TA.fix <- asreml(TA100 ~ 1 + YearPl + MaH1+and(MaH2),
#                  random = ~ ped(ID),
#                  ginverse=list(ID=Ainv),
#                  data=adata)
# 
# Cn.fix <- asreml(ES_Cn ~ 1 + PgH1+and(PgH2),
#                  random = ~ ped(ID),
#                  ginverse=list(ID=Ainv),
#                  data=adata)
# 
# M.fix <- asreml(ES_M2_lb ~ 1 + MaH1 + and(MaH2) + PgH1 + and(PgH2),
#                 random = ~ ped(ID),
#                 ginverse=list(ID=Ainv),
#                 data=adata) 

##################
# Mbest analysis #
##################
# fitted value (breeding value) = intercept(popmean) + 
#                                 ID_coefficient + 
#                                 Coef.marker(e.g., Coef.PgH1 + Coef.PgH2 + Coef.MaH1 + Coef.MaH2) +
#                                 resid
# phenotype = fitted value + residue

TA.r <- TAbest.calc(adata, TAbest, random = T)
TA.f <- TAbest.calc(adata, TA.fix, random = F)
Cn.r <- Cnbest.calc(adata, Cnbest, random = T)
Cn.f <- Cnbest.calc(adata, Cn.fix, random = F)
M.r <- Mbest.calc(adata, Mbest, random = T)
M.f <- Mbest.calc(adata, M.fix, random = F)

cor(TA.r$summary$phenotype, TA.r$summary$fit) # 0.8324748
cor(TA.f$summary$phenotype, TA.f$summary$fit) # 0.8329526
cor(Cn.r$summary$phenotype, Cn.r$summary$fit) # 0.6357701
cor(Cn.f$summary$phenotype, Cn.f$summary$fit) # 0.6131929
cor(M.r$summary$phenotype, M.r$summary$fit) # 0.8509986
cor(M.f$summary$phenotype, M.f$summary$fit) # 0.8369666

predictions <- list("TA.r" = TA.r, 
                    "TA.f" = TA.f,
                    "Cn.r" = Cn.r,
                    "Cn.f" = Cn.f,
                    "M.r" = M.r,
                    "M.f" = M.f)

save(predictions, file = "prediction.rda")


#write prediction summary to file
load(file = "prediction.rda", verbose = T)
write.func(predictions, filename = "prediction.summary.csv")
  










