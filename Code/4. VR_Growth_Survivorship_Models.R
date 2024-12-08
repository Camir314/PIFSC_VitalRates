# #Which file is the colonylevel dataframe being referenced?
# # ColTrans <- read.csv("./CSV files/ColTrans/ASRAMP23_VitalRates_colonylevel_CLEAN.csv") 
# # ColTrans <- read.csv("./CSV files/ColonyTransitions/ASRAMP23_ColonyTransitions.csv")
# loadnames=load("./Data/ColonyTransitions/Script_Step3_DataPackage.rdata");loadnames
# 
# ## Carry on ColTrans from Script 3.
# ColTrans=ColTrans %>% rename(
#   log10_SS=l10_Area_STA,
#   log10_ES=l10_Area_END,
#   log10_ESc=l10_Area_ENDann,
#   log10_SSsa=l10_SArea_STA,
#   log10_ESsa=l10_SArea_END,
#   log10_EScsa=l10_SArea_ENDann)
# 
# #Do this at SIG ModParam Level (at end of Script)
# loadHAname=load("./Data/ColonyTransitions/HA_Trans/Colony_Data_20210917_edited.rdata")
# names(ColonyLevel)
# ColTransHA=ColonyLevel %>%
#   filter(DataOrError=="DATA") %>%
#   mutate(Shape_Area_STA=StartingSize/10^4,
#          Shape_Leng_STA=StartingPerim/10^2,
#          Shape_Area_END=EndingSize/10^4,
#          Shape_Leng_END=EndingPerim/10^2,
#          l10TransitionMagnitude=log10(TransitionMagnitude),
#          StartingDate=as.character(StartingDate),
#          EndingDate=as.character(EndingDate)
#   ) %>%
#   select(ColonyID,StartingYear,EndingYear,N_t1,TransitionTypeSimple,StartingSize,EndingSize,StartingPerim,EndingPerim,
#          Island,Site,Genus_Code,StartingDate,N_t0,EndingDate,Interval,log10_SS,log10_ES,Interval_Years,log10_ESc,Mortality,Survival,
#          SiteInterval,SIG,SIS,SEC_NAME,REGION,Shape_Area_STA,Shape_Leng_STA,Shape_Area_END,Shape_Leng_END,l10TransitionMagnitude
#   ) %>%
#   rename(GENUS_CODE=Genus_Code,
#          Site_Genet=ColonyID,
#          Year_STA=StartingYear,
#          Year_END=EndingYear,
#          nPatches_END=N_t1,
#          TransitionType=TransitionTypeSimple,
#          Shape_Area_STAcm2=StartingSize,
#          Shape_Area_ENDcm2=EndingSize,
#          Shape_Leng_STAcm=StartingPerim,
#          Shape_Leng_ENDcm=EndingPerim
#   )
# 
# # setdiff(names(ColTrans),names(ColTransHA))
# # setdiff(names(ColTransHA),names(ColTrans))
# ColTransHAMAAS=full_join(ColTrans,ColTransHA)
# ColTransHAMAAS=ColTransHAMAAS %>% select(REGION,Island,SEC_NAME,Site,GENUS_CODE,Site_Genet,Interval,SIG,SIS,SiteInterval,
#                           TP_ID_STA,TP_ID_END,StartingDate,EndingDate,Year_STA,Year_END,Interval_Years,
#                           TransitionType,N_t0,nPatches_END,Mortality,Survival,
#                           Shape_Area_STA,Shape_Area_END,Shape_Leng_STA,Shape_Leng_END,SArea_STA,SArea_END,
#                           Shape_Area_STAcm2,Shape_Area_ENDcm2,Shape_Leng_STAcm,Shape_Leng_ENDcm,SArea_STAcm2,SArea_ENDcm2,
#                           log10_SS,log10_ES,log10_ESc,l10TransitionMagnitude,log10_SSsa,log10_ESsa,log10_EScsa,l10STransitionMagnitude)
# 
# table(ColTransHAMAAS$REGION)
# ColTransHAMAAS$TransitionType[ColTransHAMAAS$TransitionType=="GROWTH"]="GROW"
# write.csv(ColTransHAMAAS,"./Data/ColonyTransitions/FULLDOMAIN_19.22.23_ColonyTransitions.csv",row.names = F)


rm(list=ls())
# Loading Libraries -------------------------------------------------------
library(tidyverse)
library(lubridate)
library(stringr) 
library(pROC)
library(ggpubr)
library(sf)
library(rgeos)
library(sp)
library(patchwork)

loadnames=load("./Data/ColonyTransitions/Script_Step3_DataPackage.rdata");loadnames
#ColTrans=read.csv("./Data/ColonyTransitions/FULLDOMAIN_19.22.23_ColonyTransitions.csv")

Name="HA_MA_AS_Models"


################################################################################
######### Site + Interval Years + Genus models ###############################
########## GROWTH ! ############################################################
################################################################################
# Change this value to change the # of minimum data points at each site interval. If it is below this value then the model will not run and will skip this site/interval combination
GrowthR2_Standard <- 0.5 
MinN_Standard <- 6
PlanarGrowthTib=ColTrans %>% 
  filter(TransitionType %in% c("GROW")) %>% 
  group_by(REGION,Island,SEC_NAME,SIG,Site,Interval,GENUS_CODE,StartingDate,EndingDate,Interval_Years) %>% 
  nest() %>% 
  mutate(
    pGrowthMod=map(data,~lm(log10_ESc ~ log10_SS,data=as.data.frame(data))),
    pg.N=map(data,~length(.$log10_SS)),
    pg.int=map(pGrowthMod,~coef(.)[1]),
    pg.slp=map(pGrowthMod,~coef(.)[2]),
    pg.var=map(pGrowthMod,~(summary(.)$sigma)^2),
    pg.R2=map(pGrowthMod,~summary(.)$adj.r.squared),
    pg.AIC=map(pGrowthMod,~AIC(.)),
  ) %>% 
  unnest(c(SIG, Site,Interval,GENUS_CODE,StartingDate,Interval_Years,EndingDate,pg.N,pg.int,pg.slp,pg.var,pg.R2,pg.AIC)) %>% 
  select_if(negate(is.list)) %>% 
  mutate(pg.BadModel=ifelse(((is.na(pg.R2))|(pg.R2==1)|(pg.R2<GrowthR2_Standard)|(pg.N<MinN_Standard)),TRUE,FALSE))
PlanarGrowthDataFrame = as.data.frame(PlanarGrowthTib)

PlanarGrowthDataFrame %>% ggplot(aes(pg.N,pg.R2,size=-pg.AIC,color=pg.BadModel))+geom_point()+facet_wrap("REGION")+scale_x_sqrt()+geom_vline(xintercept = 4)

SAGrowthTib=ColTrans %>% 
  filter(TransitionType %in% c("GROW")&!is.na(log10_SSsa)) %>% 
  group_by(REGION,Island,SEC_NAME,SIG,Site,Interval,GENUS_CODE,StartingDate,EndingDate,Interval_Years) %>% 
  nest() %>% 
  mutate(
    sGrowthMod=map(data,~lm(log10_EScsa ~ log10_SSsa,data=as.data.frame(data))),
    sg.N=map(data,~length(.$log10_SSsa)),
    sg.int=map(sGrowthMod,~coef(.)[1]),
    sg.slp=map(sGrowthMod,~coef(.)[2]),
    sg.var=map(sGrowthMod,~(summary(.)$sigma)^2),
    sg.R2=map(sGrowthMod,~summary(.)$adj.r.squared),
    sg.AIC=map(sGrowthMod,~AIC(.)),
  ) %>% 
  unnest(c(SIG, Site,Interval,GENUS_CODE,StartingDate,Interval_Years,EndingDate,sg.N,sg.int,sg.slp,sg.var,sg.R2,sg.AIC)) %>% 
  select_if(negate(is.list))  %>% 
  mutate(sg.BadModel=ifelse(((is.na(sg.R2))|(sg.R2==1)|(sg.R2<GrowthR2_Standard)|(sg.N<MinN_Standard)),TRUE,FALSE))
SAGrowthDataFrame = as.data.frame(SAGrowthTib)
save(list=c("PlanarGrowthDataFrame","SAGrowthDataFrame"), file = paste0("./Data/ModelData/",Name,"_Gmodfits.rdata"))
SAGrowthDataFrame %>% ggplot(aes(sg.N,sg.R2,size=-sg.AIC,color=sg.BadModel))+geom_point()+facet_wrap("REGION")+scale_x_sqrt()


#Regional model
PlanarGrowthRegionalTib=ColTrans %>% 
  filter(TransitionType %in% c("GROW")) %>% 
  group_by(GENUS_CODE,REGION) %>% 
  nest() %>% 
  mutate(
    pGrowthMod=map(data,~lm(log10_ESc ~ log10_SS,data=as.data.frame(data))),
    pg.N=map(data,~length(.$log10_SS)),
    pg.int=map(pGrowthMod,~coef(.)[1]),
    pg.slp=map(pGrowthMod,~coef(.)[2]),
    pg.var=map(pGrowthMod,~(summary(.)$sigma)^2),
    pg.R2=map(pGrowthMod,~summary(.)$adj.r.squared),
    pg.AIC=map(pGrowthMod,~AIC(.)),
  ) %>% 
  unnest(c(GENUS_CODE,pg.N,pg.int,pg.slp,pg.var,pg.R2,pg.AIC)) %>% 
  select_if(negate(is.list)) %>% 
  mutate(pg.BadModel=ifelse(((is.na(pg.R2))|(pg.R2==1)|(pg.R2<GrowthR2_Standard)|(pg.N<MinN_Standard)),TRUE,FALSE)) 
head(PlanarGrowthRegionalTib)
PlanarGrowthRegionalDF = as.data.frame(PlanarGrowthRegionalTib)

PlanarGrowthRegionalDF %>% ggplot(aes(pg.N,pg.R2,size=-pg.AIC,color=pg.BadModel))+geom_point()+facet_wrap("REGION")+scale_x_sqrt()


SAGrowthRegionalTib=ColTrans %>% 
  filter(TransitionType %in% c("GROW")&!is.na(log10_SSsa)) %>% 
  group_by(GENUS_CODE,REGION) %>% 
  nest() %>% 
  mutate(
    sGrowthMod=map(data,~lm(log10_EScsa ~ log10_SSsa,data=as.data.frame(data))),
    sg.N=map(data,~length(.$log10_SS)),
    sg.int=map(sGrowthMod,~coef(.)[1]),
    sg.slp=map(sGrowthMod,~coef(.)[2]),
    sg.var=map(sGrowthMod,~(summary(.)$sigma)^2),
    sg.R2=map(sGrowthMod,~summary(.)$adj.r.squared),
    sg.AIC=map(sGrowthMod,~AIC(.)),
  ) %>% 
  unnest(c(GENUS_CODE,sg.N,sg.int,sg.slp,sg.var,sg.R2,sg.AIC)) %>% 
  select_if(negate(is.list))  %>% 
  mutate(sg.BadModel=ifelse(((is.na(sg.R2))|(sg.R2==1)|(sg.R2<GrowthR2_Standard)|(sg.N<MinN_Standard)),TRUE,FALSE))
head(SAGrowthRegionalTib)
SAGrowthRegionalDF = as.data.frame(SAGrowthRegionalTib)

SAGrowthRegionalDF %>% ggplot(aes(sg.N,sg.R2,size=-sg.AIC,color=sg.BadModel))+geom_point()+facet_wrap("REGION")+scale_x_sqrt()

#plot(PlanarGrowthRegionalDF$pg.R2,SAGrowthRegionalDF$sg.R2);abline(0,1)
save(list=c("PlanarGrowthRegionalDF","SAGrowthRegionalDF"), file = paste0("./Data/ModelData/",Name,"_RegionalGmodfits.rdata"))


################################################################################
######### Genera + Site + Interval Years  models ###############################
## SURVIVAL ! ##################################################################
################################################################################
#length(which(is.na(subset(ColTrans,TransitionType %in% c("GROWTH","SHRINK","MORT"))$log10_SS)))

#To annualize survival:
#Fit logistic regression at Site-Interval-Genus scale. Use this model to predict estimated survival probability for each individual.
#Annualize survival probability (raise to 1/IntervalYears) and add to ColTrans dataframe
# Drop SIG with Poor Model Fits
#Use annualized survival probability to refit logistic models for aggregate regional models 
#Run SIG models like normal using logistic regression


#Set up for pooling larger than SIG survival
wo=getOption("warn");options(warn = -1)
AnnSurvTib=ColTrans %>% 
  filter(TransitionType %in% c("GROW","MORT")&!is.na(log10_SS)) %>% 
  group_by(REGION,Island,SEC_NAME,SIG,Site,Interval,GENUS_CODE,StartingDate,EndingDate,Interval_Years) %>% 
  nest() %>% 
  mutate(
    s.N=map(data,~length(.$Survival)),
    s.ssRange=map(data,~diff(range(.$log10_SS))),
    SurvMod=map(data,~glm(Survival ~ log10_SS, family = "binomial" , data = as.data.frame(data))),
    s.int=map(SurvMod,~coef(.)[1]), 
    s.slp=map(SurvMod,~coef(.)[2]),
    s.pR2=map(SurvMod,~1-(summary(.)$deviance/summary(.)$null.deviance)),
    s.AUC=ifelse(map(data,~length(unique(.$Survival)))>1,
                 map2(data,SurvMod,~as.numeric(auc(.x$Survival,(predict(.y,newdata=select(.x,"log10_SS")))))),NA),
    s.AIC=map(SurvMod,~AIC(.)),
    SurvProbs=map(SurvMod,~predict(.,type = "response",na.action="na.pass"))#,
  ) %>% 
  ungroup() %>% 
  unnest(c(data,s.N,s.ssRange,s.int,s.slp,s.pR2,s.AUC,s.AIC,SurvProbs)) %>% 
  select_if(negate(is.list)) %>% 
  mutate(
    #Models to drop AUC == 1 or AUC < 0.50 or SIG without any mort/surv differences
    s.BadModel=ifelse((is.na(s.AUC)|(s.AUC<0.50)|(s.AUC==1)),TRUE,FALSE),
  )
AnnSurvTib$AnnSurvProbs=AnnSurvTib$SurvProbs^(1/AnnSurvTib$Interval_Years)
options(warn = wo)

#checking to make sure predict call works and we're getting actual annualized surv
# AnnSurvTib %>% filter(!is.na(s.AUC)) %>% 
#   ggplot( aes(x=s.pR2))+ 
#   geom_point(aes(y=s.AUC,size=(s.N),color=s.N>12,
#                  fill=(s.AUC<.55|s.AUC==1)),lwd=5,shape=21)+#color = "blue")+
#   # geom_point(aes(y=SurvProbs,fill=s.AUC), shape=21)+#color = "blue")+
#   # geom_point(aes(y=AnnSurvProbs,fill=s.AUC), shape=23)+#, color = "red")+
#   scale_fill_discrete(direction = -1)+
#   scale_size_area()+
#   theme_bw()
# 
# Usig=unique(AnnSurvTib$SIG)
# temp = subset(AnnSurvTib,SIG==Usig[12])
# ggplot(AnnSurvTib, aes(x=log10_SS))+ 
#   geom_point(aes(y=SurvProbs),fill="white",shape=21)+#color = "blue")+
#   geom_point(aes(y=AnnSurvProbs,fill=!s.BadModel),shape=23)+#, color = "red")+
#   facet_wrap("SIG")+theme_bw()


ColTrans=left_join(ColTrans,AnnSurvTib[,c("SIG","Site_Genet","s.N","s.AUC","s.BadModel","AnnSurvProbs")])
ColTrans=ColTrans %>% 
  mutate(TransitionMagnitudecm2=Shape_Area_ENDcm2-Shape_Area_STAcm2,
         LinearExtension.cm=TransitionMagnitudecm2/Shape_Leng_STAcm,
         AnnualLinearExtension.cm=LinearExtension.cm/Interval_Years)
save(list = "ColTrans",file="./Data/ColonyTransitions/HA_MA_AS_ColonyTransition_withAnnSurv.rdata")

ColTrans %>% filter(s.BadModel==TRUE) %>%
  ggplot(aes(x=log10_SS,y=AnnSurvProbs,fill=Site,size=s.N))+
  geom_point(shape=21,color="white")+
  facet_grid(REGION~GENUS_CODE)

### JUST REGIONAL by Genus
RegionalAnnSurvTib=ColTrans %>% 
  filter(TransitionType %in% c("GROW","MORT")&!is.na(log10_SS)) %>% 
  filter(s.BadModel==FALSE) %>% 
  group_by(GENUS_CODE,REGION) %>% 
  nest() %>% 
  mutate(
    SurvMod=map(data,~glm(AnnSurvProbs ~ log10_SS, family = "binomial" , data = as.data.frame(data))),
    NullSurvMod=map(data,~glm(AnnSurvProbs ~ 1, family = "binomial" , data = as.data.frame(data))),
    s.N=map(data,~length(.$AnnSurvProbs)),
    s.int=map(SurvMod,~coef(.)[1]), 
    s.slp=map(SurvMod,~coef(.)[2]),
    s.pR2=map(SurvMod,~1-(summary(.)$deviance/summary(.)$null.deviance)),
    s.AIC=map(SurvMod,~AIC(.)),
    s.AUC=ifelse(map(data,~length(unique(.$Survival)))>1,
                 map2(data,SurvMod,~as.numeric(auc(.x$Survival,(predict(.y,newdata=select(.x,"log10_SS")))))),NA)
  ) %>% 
  unnest(c(GENUS_CODE,s.N,s.int,s.slp, s.pR2,s.AIC,s.AUC)) %>% 
  select_if(negate(is.list)) %>% 
  mutate(
    #Models to drop AUC == 1 or AUC < 0.55 or SIG without any mort/surv differences
    s.BadModel=ifelse((is.na(s.AUC)|(s.AUC<0.55)|(s.AUC==1)),TRUE,FALSE),
  )
RegionalAnnSurvDataFrame = as.data.frame(RegionalAnnSurvTib)
head(RegionalAnnSurvDataFrame)
save(RegionalAnnSurvDataFrame, file = sprintf("./Data/ModelData/%s_RegionalSmodfits.rdata", Name))

# SIG Survival
SurvTib=ColTrans %>% 
  filter(TransitionType %in% c("GROW","MORT")&!is.na(log10_SS)) %>% 
  group_by(REGION,Island,SEC_NAME,SIG,Site,Interval,GENUS_CODE,StartingDate,EndingDate,Interval_Years) %>% 
  nest() %>% 
  mutate(
    SurvMod=map(data,~glm(Survival ~ log10_SS, family = "binomial" , data = as.data.frame(data))),
    NullSurvMod=map(data,~glm(Survival ~ 1, family = "binomial" , data = as.data.frame(data))),
    s.N=map(data,~length(.$Survival)),
    s.int=map(SurvMod,~coef(.)[1]), 
    allsurv = map(data,~length(.$Survival)),
    s.int=map(SurvMod,~coef(.)[1]),
    s.slp=map(SurvMod,~coef(.)[2]),
    s.pR2=map(SurvMod,~1-(summary(.)$deviance/summary(.)$null.deviance)),
    s.AIC=map(SurvMod,~AIC(.)),
    s.AUC=ifelse(map(data,~length(unique(.$Survival)))>1,
                 map2(data,SurvMod,~as.numeric(auc(.x$Survival,predict(.y,newdata=select(.x,"log10_SS"))))),
                 NA)
  ) %>% 
  unnest(c(SIG,Site,Interval,GENUS_CODE,StartingDate,Interval_Years,EndingDate,s.N,s.int,s.slp, s.pR2,s.AIC,s.AUC)) %>% 
  select_if(negate(is.list)) %>% 
  mutate(
    #Models to drop AUC == 1 or AUC < 0.55 or SIG without any mort/surv differences
    s.BadModel=ifelse((is.na(s.AUC)|(s.AUC<0.55)|(s.AUC==1)),TRUE,FALSE),
  )
SurvDataFrame = as.data.frame(SurvTib)
head(SurvDataFrame)
save(SurvDataFrame, file = sprintf("./Data/ModelData/%s_Smodfits.rdata", Name))

#compare regional vs SIG slope/int
GrowthPlot=ggplot()+
  geom_point(aes(pg.int,pg.slp,size=log10(pg.N),fill=REGION),alpha=.5,shape=23,data=subset(PlanarGrowthDataFrame,pg.BadModel==FALSE))+
  geom_point(aes(pg.int,pg.slp,size=log10(pg.N),fill=REGION),shape=21,data=PlanarGrowthRegionalDF)+
  facet_grid(~GENUS_CODE)
SurvPlot=ggplot()+
  geom_point(aes(s.int,s.slp,size=log10(s.N),fill=REGION),alpha=.5,shape=23,data=subset(SurvDataFrame,s.BadModel==FALSE))+
  geom_point(aes(s.int,s.slp,size=log10(s.N),fill=REGION),shape=21,data=RegionalAnnSurvDataFrame)+
  facet_grid(~GENUS_CODE)
GrowthPlot/SurvPlot

################################################################################
######### Site + Interval Years + Genus models #################################
########## INTEGRATION ! #######################################################
################################################################################
# #region lookup
# RegLU=c("NWHI","MHI","NWHI","NWHI","MHI","MHI","NWHI")
# names(RegLU)=c("FFS","HAW","KUR","LIS","MAI","OAH","PHR")
RegLU=c("NWHI","MHI","NWHI","NWHI","MHI","MHI","NWHI","MARIAN","PRIAs","MARIAN","PRIAs","MARIAN","SAMOA","MARIAN","SAMOA","MARIAN","SAMOA","SAMOA")
names(RegLU)=c("FFS","HAW","KUR","LIS","MAI","OAH","PHR","ASC","BAK", "GUA", "HOW", "MAU", "OFU", "PAG", "ROS", "SAI", "TAU", "TUT")

#SIG models
combo <- left_join(PlanarGrowthDataFrame,SurvDataFrame)
save(combo, file = sprintf("./Data/ModelData/%s_allmodfits.rdata",Name))
table(combo$pg.BadModel,combo$s.BadModel)
no_NAs <- combo %>% filter(pg.BadModel==FALSE&s.BadModel==FALSE)

GrowthSurv_SIG <- data.frame()

GrowthSurv_SIG <- no_NAs[order(no_NAs$Site),]
GrowthSurv_SIG$StartingYear=as.numeric(year(GrowthSurv_SIG$StartingDate))
#GrowthSurv_SIG$SEC_NAME#=Site2Sec_ONLY[match(GrowthSurv_SIG$Site,Site2Sec_ONLY$Site),"SEC_NAME"]

#sort(table(GrowthSurv_SIG$Site,GrowthSurv_SIG$GENUS_CODE,GrowthSurv_SIG$StartingYear))

#add SfM site-level stock recruitment data
GrowthSurvRec_SIG = left_join(GrowthSurv_SIG,
                              RecSFMDataFrame[,c("Site","GENUS_CODE","EndingDate","Nrec","RecSFM_p_SiteAdcm2_Yr")],
                              by=c("Site","GENUS_CODE","StartingDate"="EndingDate"))
#View(GrowthSurvRec_SIG[,c("SIG","StartingDate","RecSFM_p_SiteAdcm2_Yr")])

#add REA sector-level stock recruitment data
#Tweak ANALYSIS YEAR values to match recruit data
GrowthSurvRec_SIG$ANALYSIS_YEAR=as.character(GrowthSurvRec_SIG$StartingYear)
GrowthSurvRec_SIG$ANALYSIS_YEAR[which(GrowthSurvRec_SIG$ANALYSIS_YEAR=="2013")]="2013-15"
GrowthSurvRec_SIG$ANALYSIS_YEAR[which(GrowthSurvRec_SIG$ANALYSIS_YEAR=="2015")]="2013-15"
GrowthSurvRec_SIG$ANALYSIS_YEAR[which(GrowthSurvRec_SIG$ANALYSIS_YEAR%in%c("2017","2018")&GrowthSurvRec_SIG$SEC_NAME%in%c("Baker","Howland"))]="2017-18"

# GrowthSurvRec_SIG %>% arrange(SEC_NAME,StartingYear) %>%  select(c("SEC_NAME","GENUS_CODE","ANALYSIS_YEAR")) %>%
#   left_join(Jsec_P_Asec[,c("ANALYSIS_SEC","ANALYSIS_YEAR","GENUS_CODE","MN.RecVal")],by=c("SEC_NAME"="ANALYSIS_SEC","GENUS_CODE","ANALYSIS_YEAR"))
# Jsec_P_Asec %>% ungroup() %>%  select(c("ANALYSIS_SEC","ANALYSIS_YEAR")) %>% distinct()

ModelParams_SIG= left_join(GrowthSurvRec_SIG,
                           Jsec_P_Asec[,c("ANALYSIS_SEC","GENUS_CODE","ANALYSIS_YEAR","N.RecVal","MN.RecVal","LOW_CI95.RecVal","HIGH_CI95.RecVal")],
                           by=c("SEC_NAME"="ANALYSIS_SEC","GENUS_CODE","ANALYSIS_YEAR"))

ModelParams_SIG[,c("SIG","SEC_NAME","MN.RecVal")] %>% filter(is.na(MN.RecVal))

ModelParams_SIG= ModelParams_SIG %>% rename(
  N.Rec_Site=Nrec,
  MN.RecVal_Site=RecSFM_p_SiteAdcm2_Yr,
  N.RecVal_Sec=N.RecVal,
  MN.RecVal_Sec=MN.RecVal,
  LOW_CI95.RecVal_Sec=LOW_CI95.RecVal,
  HIGH_CI95.RecVal_Sec=HIGH_CI95.RecVal
)


#add SSSS (Scleractinia) sector rec values
Rec_SSSS=subset(Jsec_P_Asec,GENUS_CODE=="SSSS")
names(Rec_SSSS)[21:26]=paste0(names(Rec_SSSS)[21:26],"_Sec_SSSS")
ModelParams_SIG= left_join(ModelParams_SIG,
                           Rec_SSSS[,c("ANALYSIS_SEC","ANALYSIS_YEAR",
                                       "N.RecVal_Sec_SSSS","MN.RecVal_Sec_SSSS","LOW_CI95.RecVal_Sec_SSSS","HIGH_CI95.RecVal_Sec_SSSS")],
                           by=c("SEC_NAME"="ANALYSIS_SEC","ANALYSIS_YEAR"))
# ModelParams_SIG$ISLAND=substr(ModelParams_SIG$Site,5,7)
# ModelParams_SIG$REGION=RegLU[ModelParams_SIG$ISLAND]
ModelParams_SIG= left_join(ModelParams_SIG,
                           RecVal_Sec_Dists[,c("REGION","GENUS_CODE",
                                               "MD.RecVal_Sec_All","MD.CI95_LO.RecVal_Sec_All","MD.CI95_HI.RecVal_Sec_All")],
                           by=c("REGION","GENUS_CODE"))

table(!is.na(ModelParams_SIG$MN.RecVal_Site))
table(!is.na(ModelParams_SIG$MN.RecVal_Sec))
table(!is.na(ModelParams_SIG$MN.RecVal_Site)&!is.na(ModelParams_SIG$MN.RecVal_Sec))
table(!is.na(ModelParams_SIG$MD.RecVal_Sec_All))

# write.csv(RecSFMDataFrame,file = "./Data/ModelData/HAMAAS_Site_Recruits_SFM.csv")
# write.csv(Jsec_P_Asec,file = "./Data/ModelData/HAMAAS_Sector_Recruits_REA.csv")
# write.csv(ModelParams_SIG,file = "./Data/ModelData/HAMAAS_ModelParams_SIG.csv")



ModelParams_SIG=
  ModelParams_SIG %>%
    mutate(SecRecVal=coalesce(MN.RecVal_Sec,MN.RecVal_Sec_SSSS,MD.RecVal_Sec_All),
           SecRecValSourceS=ifelse(SecRecVal==MN.RecVal_Sec,"SEC","NA"),
           SecRecValSourceSS=ifelse(SecRecVal==MN.RecVal_Sec_SSSS,"SEC_SSSS","NA"),
           SecRecValSourceALL=ifelse(SecRecVal==MD.RecVal_Sec_All,"SEC_REG","NA"),
           SecRecValSource=coalesce(SecRecValSourceS,SecRecValSourceSS,SecRecValSourceALL)
    ) %>% select(!c(SecRecValSourceS,SecRecValSourceSS,SecRecValSourceALL))
      
save(ModelParams_SIG, file = sprintf("./Data/ModelData/%s_allmodfits_noNAs.rdata",Name))

ModelParams_SIG %>% 
  group_by(Site,Interval) %>% 
  summarize(N=sum(pg.N)) %>% 
  pivot_wider(names_from = Interval, values_from = N) %>% 
  print(n=nrow(.))



################################################################################
######### Regional + Genus models #################################
########## INTEGRATION ! #######################################################
################################################################################

#Regional models (genus and region)
combo_regional <- left_join(PlanarGrowthRegionalDF,RegionalAnnSurvDataFrame)#,SAGrowthRegionalDF)
combo_regional <- (combo_regional)%>% filter(pg.BadModel==FALSE&s.BadModel==FALSE)

#add SfM site-level stock recruitment data
ModelParams_Regional = left_join(combo_regional,RecSFMSummary, by= c("GENUS_CODE","REGION") )

#add REA sector-level stock recruitment data
ModelParams_Regional = left_join(ModelParams_Regional,Regional_SectorStockRec,by=c("GENUS_CODE"="GENUS_CODE","REGION"))

save(ModelParams_Regional, file = sprintf("./Data/ModelData/%s_RegionalModFits_noBADs.rdata",Name))



