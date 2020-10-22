library(lattice)
library(asreml)
source(file = "S1a_source.R")
###############################################################################
# select prediction models for acidity (TA), crispness (Cn) and firmness (M2)
###############################################################################
# Significance test for a random variable
# # comparison between M02 and M01 (< 0.1 -> significant. Diff = sig at 0.05 if < 0.1)
# 1-pchisq(2*(M02$logl-M01$logl),(M02.df-M01.df))
########################################
#### Step 1. read & formating data  ####
########################################
# load data
load(file = "data.frmt.rda")
adata <- input$data_frmt
Aniv <- input$Ainv
ped_in <- input$ped_in

####################################
#####  1. test for acidity   #######
####################################
# TA01-- simplest model.
TA01 <- asreml(TA10 ~ 1,
              random = ~ ped(ID),
              ginverse=list(ID=Ainv),
              data=adata)

print(log.TA01 <- summary(TA01)$logl) #-117.2106
print(var.TA01 <- summary(TA01)$var)
# gamma component std.error  z.ratio constraint
# gamma component std.error  z.ratio constraint
# ped(ID)!ped 1.379992  58.73872  27.34256 2.148252   Positive
# R!variance  1.000000  42.56455  13.93137 3.055304   Positive
print(VA <- var.TA01$component[1]) # 0.5874
print(VE <- var.TA01$component[2]) # 0.4256
print(h2 <- VA/(VA+VE)) # 0.58


######################################
# TA02-- TA10 ~ 1 + YearPl
######################################
# Better than TA01, YearPl significant
TA02 <- asreml(TA10 ~ 1 + YearPl,
              random = ~ ped(ID),
              ginverse=list(ID=Ainv),
              data=adata)

# wald test - for fixed effect
wald.asreml(TA02)
# Wald tests for fixed effects
# 
# Response: TA10
# 
# Terms added sequentially; adjusted for those above
# 
# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1    74.593        126.906 < 2.2e-16 ***
#   YearPl         1     4.994          8.496  0.003559 ** 
#   residual (MS)        0.588                             
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

print(log.TA02 <- summary(TA02)$logl) #-116.3411

#################################
# TA03 -- MaH1 MaH2 random effect
#################################
# Better than TA02, MaH1 MaH2 significant
TA03 <- asreml(TA10 ~ 1 + YearPl,
              random = ~ ped(ID)+MaH1+and(MaH2),
              ginverse=list(ID=Ainv),
              data=adata)
print(log.TA03 <- summary(TA03)$logl) #-105.4084
print(sig.TA03 <- sig(TA03, TA02, 1)) 
# sig     logl.diff      m.t.logl    m.ori.logl 
# 2.924749e-06  1.093265e+01 -1.054084e+02 -1.163411e+02

############################################################
# TA04 -- PgH1 PgH2 random effect
############################################################
### not better than TA03, PgH1, PgH2 not significant
TA04 <- asreml(TA10 ~ 1 + YearPl,
              random = ~ ped(ID)+MaH1+and(MaH2)+PgH1 + and(PgH2),
              ginverse=list(ID=Ainv),
              data=adata)
print(log.TA04 <- summary(TA04)$logl) #-105.4084
print(sig.TA04 <- sig(TA04, TA03, 1)) 
# sig     logl.diff      m.t.logl    m.ori.logl 
# 1.000000e+00 -7.003306e-06 -8.399331e+02 -8.399331e+02 

############################################################
# TA05 -- TA100~ YearPl, random = ped(ID) + MaH1+and(MaH2)+ ACS1H1 + and(ACS1H2)
############################################################
# not better than M03, ACS1 not significant
TA05 <- asreml(TA10 ~ 1 + YearPl,
              random = ~ ped(ID) + MaH1+and(MaH2)+ ACS1H1 + and(ACS1H2),
              ginverse=list(ID=Ainv),
              data=adata)

print(log.TA05 <- summary(TA05)$logl) #-105.2118
print(sig.TA05 <- sig(TA05, TA03, 1)) 
# sig    logl.diff     m.t.logl   m.ori.logl 
# 0.5306055    0.1966159 -105.2118138 -105.4084297

############################################################
# TA06 -- TA10~ YearPl, random = ped(ID) + MaH1+and(MaH2)+ ACO1H1 + and(ACO1H2)
############################################################
# not better than TA03, ACO1 not significant
TA06 <- asreml(TA10 ~ 1 + YearPl,
              random = ~ ped(ID) + MaH1+and(MaH2)+ ACO1H1 + and(ACO1H2),
              ginverse=list(ID=Ainv),
              data=adata)
print(log.TA06 <- summary(TA06)$logl) #-105.4085
print(sig.TA06 <- sig(TA06, TA03, 1)) 
# sig     logl.diff      m.t.logl    m.ori.logl 
# 1.000000e+00 -9.099238e-05 -1.054085e+02 -1.054084e+02

####################################
#####  2. test for crispness #######
####################################
# Cn01 -- simplest model
Cn01 <- asreml(ES_Cn ~ 1,
              random = ~ ped(ID),
              ginverse=list(ID=Ainv),
              data=adata)

log.Cn01 <- summary(Cn01)$logl #-1487.822
var.Cn01 <- summary(Cn01)$var
#             gamma component std.error  z.ratio constraint
# ped(ID)!ped 0.3874644  1236.986  725.0469 1.706077   Positive
# R!variance  1.0000000  3192.514  465.1616 6.863236   Positive
VA <- var.Cn01$component[1] # 1236.986
VE <- var.Cn01$component[2] # 3192.514
h2 <- VA/(VA+VE) # 0.2792608

#############################################
# Cn02 -- YearPl
#############################################
# Cn01 is better, YearPl not significant
Cn02 <- asreml(ES_Cn ~ 1 + YearPl,
              random = ~ ped(ID),
              ginverse=list(ID=Ainv),
              data=adata)

# wald test
wald.asreml(Cn02)
# Wald tests for fixed effects
# Response: ES_Cn
# Terms added sequentially; adjusted for those above
# 
# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1    162246         51.671 6.564e-13 ***
#   YearPl         1        27          0.009     0.926    
# residual (MS)         3140                             
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#############################
# Cn04 -- Ma allele effect
##############################
# Cn01 is better, Ma is not significant
Cn04 <- asreml(ES_Cn ~ 1,
              random = ~ ped(ID)+MaH1+and(MaH2),
              ginverse=list(ID=Ainv),
              data=adata)
print(sig.Cn04 <- sig(Cn04, Cn01, 1))
# sig     logl.diff      m.t.logl    m.ori.logl 
# 0.2111021     0.7819260 -1487.0402937 -1487.8222198 

######################################
# Cn05 -- ACS allele effect
######################################
# Cn01 is better, ACS not significant
Cn05 <- asreml(ES_Cn ~ 1,
              random = ~ ped(ID) + ACS1H1 + and(ACS1H2),
              ginverse=list(ID=Ainv),
              data=adata)
print(sig.Cn05 <- sig(Cn05, Cn01, 1))
# sig     logl.diff      m.t.logl    m.ori.logl 
# 0.6383220     0.1104717 -1487.7117480 -1487.8222198

############################################################
# Cn06 -- ACO allele effect
# Cn01 is better. ACO not significant
Cn06 <- asreml(ES_Cn ~ 1,
              random = ~ ped(ID) + ACO1H1 + and(ACO1H2),
              ginverse=list(ID=Ainv),
              data=adata)

print(sig.Cn06 <- sig(Cn06, Cn01, 1))
# sig     logl.diff      m.t.logl    m.ori.logl 
# 0.4449180     0.2917829 -1487.5304369 -1487.8222198

#############################################
# Cn08 -- Pg randome effect
# Cn08 better than Cn01, Pg is significant
Cn08 <- asreml(ES_Cn ~ 1,
              random = ~ ped(ID) + PgH1+and(PgH2),
              ginverse=list(ID=Ainv),
              data=adata)

print(sig.Cn08 <- sig(Cn08, Cn01, 1))
# sig     logl.diff      m.t.logl    m.ori.logl 
# 0.01428625  3.001211e+00 -1.484821e+03 -1.487822e+03


####################################
#####  3. test for firmness  #######
####################################
# M01-- simplest model
M01 <- asreml(ES_M2_lb ~ 1,
              random = ~ ped(ID),
              ginverse=list(ID=Ainv),
              data=adata)
print(log.M01 <- summary(M01)$logl) #-587.6946
print(var.M01 <- summary(M01)$var)
# gamma component std.error  z.ratio constraint
# ped(ID)!ped 1.413768 11.275967  5.282130 2.134738   Positive
# R!variance  1.000000  7.975827  2.684592 2.970964   Positive
print(VA <- var.M01$component[1]) # 11.27597
print(VE <- var.M01$component[2]) # 7.975827
print(h2 <- VA/(VA+VE)) # 0.5857099

#############################
# M02-- simplest model+YearPl
# no bettter than M01, YearPl not significant
M02 <- asreml(ES_M2_lb ~ 1 + YearPl,
              random = ~ ped(ID),
              ginverse=list(ID=Ainv),
              data=adata)
print(log.M02 <- summary(M02)$logl) #-586.952
# wald test
# wald.asreml(M02)
# Terms added sequentially; adjusted for those above
# 
# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1   1957.00        201.370    <2e-16 ***
#   YearPl         1     21.54          2.216    0.1366    
# residual (MS)         9.72                             
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#########################
# M04 -- Ma allele effect
# better than M03, Ma allele effect significant
M04 <- asreml(ES_M2_lb ~ 1,
              random = ~ ped(ID) + MaH1 + and(MaH2),
              ginverse=list(ID=Ainv),
              data=adata)

print(log.M04 <- summary(M04)$logl) #-586.2402
print(sig.M04 <- sig(M04, M01, 1)) 
# sig     logl.diff      m.t.logl    m.ori.logl 
# 0.08809423    1.45442907 -586.24019741 -587.69462648

#########################
# M05 -- Pg allele effect
# best model so far, Pg effect significant
M05 <- asreml(ES_M2_lb ~ 1,
              random = ~ ped(ID) + MaH1 + and(MaH2) + PgH1 + and(PgH2),
              ginverse=list(ID=Ainv),
              data=adata)

print(log.M05 <- summary(M05)$logl) #-580.4948
print(sig.M05 <- sig(M05, M04, 1)) 
# sig     logl.diff      m.t.logl    m.ori.logl 
# 6.994320e-04  5.745378e+00 -5.804948e+02 -5.862402e+02

##########################
# M06 -- ACS allele effect
# no better than M05, ACS not significant
M06 <- asreml(ES_M2_lb ~ 1,
              random = ~ ped(ID) + MaH1 + and(MaH2) + PgH1 + and(PgH2) + ACS1H1 + and(ACS1H2),
              ginverse=list(ID=Ainv),
              data=adata)

print(log.M06 <- summary(M06)$logl) #-580.3388
print(sig.M06 <- sig(M06, M05, 1)) 
# sig    logl.diff     m.t.logl   m.ori.logl 
# 0.5764304    0.1560205 -580.3387991 -580.4948196

############################
# M07 -- ACO allele effect
# no better than M05, ACO not significant
M07 <- asreml(ES_M2_lb ~ 1,
              random = ~ ped(ID) + MaH1 + and(MaH2) + PgH1 + and(PgH2) + ACO1H1 + and(ACO1H2),
              ginverse=list(ID=Ainv),
              data=adata)
print(log.M07 <- summary(M07)$logl) #-580.4948
print(sig.M07 <- sig(M07, M05, 1)) 
# sig     logl.diff      m.t.logl    m.ori.logl 
# 1.000000e+00 -5.321444e-06 -5.804948e+02 -5.804948e+02 



###########################################
#######       best models        ##########
###########################################
TAbest <- TA03 <- asreml(TA10 ~ 1 + YearPl,
                         random = ~ ped(ID)+MaH1+and(MaH2),
                         ginverse=list(ID=Ainv),
                         data=adata)

Cnbest <- Cn08 <- asreml(ES_Cn ~ 1,
               random = ~ ped(ID) + PgH1+and(PgH2),
               ginverse=list(ID=Ainv),
               data=adata)

Mbest <- M05 <- asreml(ES_M2_lb ~ 1,
              random = ~ ped(ID) + MaH1 + and(MaH2) + PgH1 + and(PgH2),
              ginverse=list(ID=Ainv),
              data=adata) 


# fixed effect models
TA.fix <- asreml(TA10 ~ 1 + YearPl + MaH1+and(MaH2),
                 random = ~ ped(ID),
                 ginverse=list(ID=Ainv),
                 data=adata)
wald.asreml(TA.fix)

# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1    45.631         96.358 < 2.2e-16 ***
#   YearPl         1     2.750          5.808   0.01595 *  
#   MaH1           6    17.979         37.966 1.141e-06 ***
#   residual (MS)        0.474                             
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# insignificant Pg
# TA.fix.temp <- asreml(TA10 ~ 1 + YearPl + PgH1+and(PgH2),
#                     random = ~ ped(ID),
#                     ginverse=list(ID=Ainv),
#                     data=adata)
# wald.asreml(TA.fix.temp)

# insignificant ACS
# TA.fix.temp <- asreml(TA10 ~ 1 + YearPl + ACS1H1+and(ACS1H2),
#                     random = ~ ped(ID),
#                     ginverse=list(ID=Ainv),
#                     data=adata)
# wald.asreml(TA.fix.temp)

# insignificant ACO
# TA.fix.temp <- asreml(TA10 ~ 1 + YearPl + ACO1H1+and(ACO1H2),
#                       random = ~ ped(ID),
#                       ginverse=list(ID=Ainv),
#                       data=adata)
# wald.asreml(TA.fix.temp)


Cn.fix <- asreml(ES_Cn ~ 1 + PgH1+and(PgH2),
                 random = ~ ped(ID),
                 ginverse=list(ID=Ainv),
                 data=adata)
wald.asreml(Cn.fix)
# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1    280984         84.363   < 2e-16 ***
#   PgH1           2     41809         12.553   0.00188 ** 
#   residual (MS)         3331                             
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

M.fix <- asreml(ES_M2_lb ~ 1 + MaH1 + and(MaH2) + PgH1 + and(PgH2),
                random = ~ ped(ID),
                ginverse=list(ID=Ainv),
                data=adata) 
wald.asreml(M.fix)
# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1   2005.18        221.975 < 2.2e-16 ***
#   MaH1           6    159.27         17.631  0.007223 ** 
#   PgH1           2    168.26         18.626 9.023e-05 ***
#   residual (MS)         9.03                             
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


Model.select <- list("TAbest" = TAbest, "Cnbest" = Cnbest, "Mbest" = Mbest, "TA.fix" = TA.fix, "Cn.fix" = Cn.fix, "M.fix" = M.fix)
save(Model.select, file = "Model.select.rda")
