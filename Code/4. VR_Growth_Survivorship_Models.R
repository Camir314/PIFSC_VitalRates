#Which file is the colonylevel dataframe being referenced?
ColonyLevel <- read.csv("./CSV files/ColonyLevel/ASRAMP23_VitalRates_colonylevel_CLEAN.csv") 
ColonyLevel <- read.csv("./CSV files/ColonyTransitions/ASRAMP23_ColonyTransitions.csv")

  ################################################################################
  ######### Site + Interval Years + Genus models ###############################
  ########## GROWTH ! ############################################################
  ################################################################################
  # Change this value to change the # of minimum data points at each site interval. If it is below this value then the model will not run and will skip this site/interval combination
  #minTrans <- 20 
  
  GrowthTib=ColonyLevel %>% 
    filter(TransitionTypeSimple %in% c("GROWTH","SHRINK")) %>% 
    group_by(SIG,Site,Interval,Genus_Code,StartingDate,EndingDate,Interval_Years) %>% 
    nest() %>% 
    mutate(
      GrowthMod=map(data,~lm(log10_ESc ~ log10_SS,data=as.data.frame(data))),
      g.N=map(data,~length(.$log10_SS)),
      g.int=map(GrowthMod,~coef(.)[1]),
      g.slp=map(GrowthMod,~coef(.)[2]),
      g.var=map(GrowthMod,~(summary(.)$sigma)^2),
      g.R2=map(GrowthMod,~summary(.)$adj.r.squared),
      g.AIC=map(GrowthMod,~AIC(.)),
    ) %>% 
    unnest(c(SIG, Site,Interval,Genus_Code,StartingDate,Interval_Years,EndingDate,g.N,g.int,g.slp,g.var,g.R2,g.AIC)) %>% 
    select_if(negate(is.list)) 
  
  GrowthDataFrame = as.data.frame(GrowthTib)
  #save(GrowthDataFrame, file = paste0("data/",Name,"_Gmodfits.rdata"))

  
  #Regional model
  GrowthRegionalTib=ColonyLevel %>% 
    filter(TransitionTypeSimple %in% c("GROWTH","SHRINK")) %>% 
    group_by(Genus_Code,REGION) %>% 
    nest() %>% 
    mutate(
      GrowthMod=map(data,~lm(log10_ESc ~ log10_SS,data=as.data.frame(data))),
      g.N=map(data,~length(.$log10_SS)),
      g.int=map(GrowthMod,~coef(.)[1]),
      g.slp=map(GrowthMod,~coef(.)[2]),
      g.var=map(GrowthMod,~(summary(.)$sigma)^2),
      g.R2=map(GrowthMod,~summary(.)$adj.r.squared),
      g.AIC=map(GrowthMod,~AIC(.)),
    ) %>% 
    unnest(c(Genus_Code,g.N,g.int,g.slp,g.var,g.R2,g.AIC)) %>% 
    select_if(negate(is.list)) 
  head(GrowthRegionalTib)
  GrowthRegionalDF = as.data.frame(GrowthRegionalTib)
  #save(GrowthRegionalDF, file = paste0("data/",Name,"_RegionalGmodfits.rdata"))
  
  ColonyLevel$Circularity=(4*pi*ColonyLevel$StartingSize)/(ColonyLevel$StartingPerim^2)
  ColonyLevel$NfragBin=cut(ColonyLevel$N_t0,breaks=c(0,1,3,10,100)+.5,labels=c("1","2-3","4-9","10+"))
  ColonyLevel$CircularityBin=cut(ColonyLevel$Circularity,breaks=seq(0,1,length.out = 5))
  ColonyLevel$DailyPropRate_E=ColonyLevel$AnnualPropRate_E^(1/365)
  ColonyLevel$MonthlyPropRate_E=ColonyLevel$AnnualPropRate_E^(1/12)
  
  ColonyLevel %>% 
    filter(TransitionTypeSimple %in% c("GROWTH")) %>% 
    filter(Genus_Code=="POSP"&REGION=="MHI") %>%
    ggplot(aes(x=StartingSummedPatchDiam,
               y=100*(MonthlyPropRate_E-1),color=NfragBin))+#y=log10(TransitionRate_L)
    geom_point(alpha=.1)+
    stat_smooth(method="loess",span=1.2)+
    scale_color_brewer(name="No. Patches\n in Colony",palette = "Set1")+
    ylab("Monthly Percentage Growth (%)")+
    scale_x_log10()+
    #scale_y_continuous(breaks=-3:3,labels=c(.001,.01,.1,1,10,100,1000))+
    #scale_y_continuous(limits=c(0,10))+
    #geom_abline()+
    geom_hline(yintercept = 0)+
    # facet_wrap("NfragBin")+
    theme_bw()
  
  PatchLevel %>% 
    filter(TransitionType %in% c("GROWTH","SHRINK")) %>% 
    filter(Genus_Code=="POSP") %>%
    ggplot(aes(x=StartingSize,y=log2(1+PercentChange)))+#y=log10(TransitionRate_L)
    geom_point(alpha=.1)+
    stat_smooth(method="loess",span=1.2)+
    scale_color_brewer(palette = "Set1")+
    #scale_y_continuous(breaks=-3:3,labels=c(.001,.01,.1,1,10,100,1000))+
    #scale_y_log10(limits=c())+
    scale_x_log10()+
    #geom_abline()+
    geom_hline(yintercept = 1)+
    # facet_wrap("NfragBin")+
    theme_bw()
  
  
  ################################################################################
  ######### Genera + Site + Interval Years  models ###############################
  ## SURVIVAL ! ##################################################################
  ################################################################################
  #length(which(is.na(subset(ColonyLevel,TransitionTypeSimple %in% c("GROWTH","SHRINK","MORT"))$log10_SS)))
  
  #To annualize survival:
  #Fit logistic regression at Site-Interval-Genus scale. Use this model to predict estimated survival probability for each individual.
  #Annualize survival probability (raise to 1/IntervalYears) and add to ColonyLevel dataframe
  #Use annualized survival probability to refit logistic models for aggregate regional models 
  #Run SIG models like normal using logistic regression
  
  
  #Set up for pooling larger than SIG survival
  wo=getOption("warn");options(warn = -1)
  AnnSurvTib=ColonyLevel %>% 
    filter(TransitionTypeSimple %in% c("GROWTH","SHRINK","MORT")&!is.na(log10_SS)) %>% 
    group_by(SIG,Site,Interval,Genus_Code,StartingDate,EndingDate,Interval_Years) %>% 
    nest() %>% 
    mutate(
      s.N=map(data,~length(.$Survival)),
      SurvMod=map(data,~glm(Survival ~ log10_SS, family = "binomial" , data = as.data.frame(data))),
      SurvProbs=map(SurvMod,~predict(.,type = "response",na.action="na.pass"))#,
    ) %>% 
    ungroup() %>% 
    unnest(c(data,s.N,SurvProbs)) %>% 
    select_if(negate(is.list)) 
  AnnSurvTib$AnnSurvProbs=AnnSurvTib$SurvProbs^(1/AnnSurvTib$Interval_Years)
  options(warn = wo)
  
  #checking to make sure predict call works and we're getting actual annualized surv
  Usig=unique(AnnSurvTib$SIG)
  temp = subset(AnnSurvTib,SIG==Usig[12])
  ggplot(temp, aes(x=log10_SS))+ geom_point(aes(y=SurvProbs), color = "blue")+ geom_point(aes(y=AnnSurvProbs), color = "red")
  
  
  ColonyLevel_ap=left_join(ColonyLevel,AnnSurvTib[,c("SIG","ColonyID","s.N","AnnSurvProbs")])
  
  ggplot(ColonyLevel_ap,aes(x=log10_SS,y=AnnSurvProbs,fill=Site,size=s.N))+
    geom_point(shape=21,color="white")+
    facet_grid(Island~Genus_Code)
  
  ### JUST REGIONAL by Genus
  RegionalAnnSurvTib=ColonyLevel_ap %>% 
    filter(TransitionTypeSimple %in% c("GROWTH","SHRINK","MORT")) %>% 
    group_by(Genus_Code,REGION) %>% 
    nest() %>% 
    mutate(
      SurvMod=map(data,~glm(AnnSurvProbs ~ log10_SS, family = "binomial" , data = as.data.frame(data))),
      NullSurvMod=map(data,~glm(AnnSurvProbs ~ 1, family = "binomial" , data = as.data.frame(data))),
      s.N=map(data,~length(.$AnnSurvProbs)),
      s.int=map(SurvMod,~coef(.)[1]), 
      allsurv = map(data,~length(.$AnnSurvProbs)),
      s.int=map(SurvMod,~coef(.)[1]),
      s.slp=map(SurvMod,~coef(.)[2]),
      s.pR2=map(SurvMod,~1-(summary(.)$deviance/summary(.)$null.deviance)),
      s.AIC=map(SurvMod,~AIC(.))
    ) %>% 
    unnest(c(Genus_Code,s.N,s.int,s.slp, s.pR2,s.AIC)) %>% 
    select_if(negate(is.list)) 
  RegionalAnnSurvDataFrame = as.data.frame(RegionalAnnSurvTib)
  head(RegionalAnnSurvDataFrame)
  #save(RegionalAnnSurvDataFrame, file = sprintf("data/%s_RegionalSmodfits.rdata", Name))
  
  # SIG Survival
  SurvTib=ColonyLevel %>% 
    filter(TransitionTypeSimple %in% c("GROWTH","SHRINK","MORT")) %>% 
    group_by(SIG,Site,Interval,Genus_Code,StartingDate,EndingDate,Interval_Years) %>% 
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
      s.AIC=map(SurvMod,~AIC(.))
    ) %>% 
    unnest(c(SIG,Site,Interval,Genus_Code,StartingDate,Interval_Years,EndingDate,s.N,s.int,s.slp, s.pR2,s.AIC)) %>% 
    select_if(negate(is.list)) 
  SurvDataFrame = as.data.frame(SurvTib)
  head(SurvDataFrame)
  #save(SurvDataFrame, file = sprintf("data/%s_Smodfits.rdata", Name))
  
  #compare regional vs SIG slope/int
  ggplot()+
    geom_point(aes(s.int,s.slp,size=s.N,shape=Genus_Code),color="red",data=RegionalAnnSurvDataFrame)+
    geom_point(aes(s.int,s.slp,size=s.N,shape=Genus_Code),data=subset(SurvDataFrame,s.N>=0))+
    facet_wrap("Genus_Code")
  
  ################################################################################
  ######### Site + Interval Years + Genus models #################################
  ########## INTEGRATION ! #######################################################
  ################################################################################
  #region lookup
  RegLU=c("NWHI","MHI","NWHI","NWHI","MHI","MHI","NWHI")
  names(RegLU)=c("FFS","HAW","KUR","LIS","MAI","OAH","PHR")
  
  #SIG models
  combo <- left_join(GrowthDataFrame,SurvDataFrame)
  #save(combo, file = sprintf("data/%s_allmodfits.rdata",Name))
  no_NAs <- na.omit(combo)
  GrowthSurv_SIG <- data.frame()
  
  GrowthSurv_SIG <- no_NAs[order(no_NAs$Site),]
  GrowthSurv_SIG$StartingYear=as.numeric(year(GrowthSurv_SIG$StartingDate))
  GrowthSurv_SIG$SEC_NAME=Site2Sec_ONLY[match(GrowthSurv_SIG$Site,Site2Sec_ONLY$Site),"SEC_NAME"]
  
  #sort(table(GrowthSurv_SIG$Site,GrowthSurv_SIG$Genus_Code,GrowthSurv_SIG$StartingYear))
  
  #add SfM site-level stock recruitment data
  GrowthSurvRec_SIG = left_join(GrowthSurv_SIG,
                                RecSFMDataFrame[,c("Site","Genus_Code","EndingDate","Nrec","RecSFM_p_SiteAdcm2_Yr")],
                                by=c("Site","Genus_Code","StartingDate"="EndingDate"))
  #add REA sector-level stock recruitment data
  ModelParams_SIG= left_join(GrowthSurvRec_SIG,
                             Jsec_P_Asec[,c("ANALYSIS_SEC","GENUS_CODE"," OBS_YEAR","N.RecVal","MN.RecVal","LOW_CI95.RecVal","HIGH_CI95.RecVal")],
                             by=c("SEC_NAME"="ANALYSIS_SEC","Genus_Code"="GENUS_CODE","StartingYear"="OBS_YEAR"))
  names(ModelParams_SIG)[which(names(ModelParams_SIG)=="Nrec")]="N.Rec_Site"
  names(ModelParams_SIG)[which(names(ModelParams_SIG)=="RecSFM_p_SiteAdcm2_Yr")]="MN.RecVal_Site"
  names(ModelParams_SIG)[which(names(ModelParams_SIG)=="N.RecVal")]="N.RecVal_Sec"
  names(ModelParams_SIG)[which(names(ModelParams_SIG)=="MN.RecVal")]="MN.RecVal_Sec"
  names(ModelParams_SIG)[which(names(ModelParams_SIG)=="LOW_CI95.RecVal")]="LOW_CI95.RecVal_Sec"
  names(ModelParams_SIG)[which(names(ModelParams_SIG)=="HIGH_CI95.RecVal")]="HIGH_CI95.RecVal_Sec"
  #add SSSS (Scleractinia) sector rec values
  Rec_SSSS=subset(Jsec_P_Asec,GENUS_CODE=="SSSS")
  names(Rec_SSSS)[21:26]=paste0(names(Rec_SSSS)[21:26],"_Sec_SSSS")
  ModelParams_SIG= left_join(ModelParams_SIG,
                             Rec_SSSS[,c("ANALYSIS_SEC","OBS_YEAR",
                                         "N.RecVal_Sec_SSSS","MN.RecVal_Sec_SSSS","LOW_CI95.RecVal_Sec_SSSS","HIGH_CI95.RecVal_Sec_SSSS")],
                             by=c("SEC_NAME"="ANALYSIS_SEC","StartingYear"="OBS_YEAR"))
  ModelParams_SIG$ISLAND=substr(ModelParams_SIG$Site,1,3)
  ModelParams_SIG$REGION=RegLU[ModelParams_SIG$ISLAND]
  ModelParams_SIG= left_join(ModelParams_SIG,
                             RecVal_Sec_Dists[,c("REGION","GENUS_CODE",
                                                 "MD.RecVal_Sec_All","MD.CI95_LO.RecVal_Sec_All","MD.CI95_HI.RecVal_Sec_All")],
                             by=c("REGION"="REGION","Genus_Code"="GENUS_CODE"))
  
  table(!is.na(ModelParams_SIG$MN.RecVal_Site))
  table(!is.na(ModelParams_SIG$MN.RecVal_Sec))
  table(!is.na(ModelParams_SIG$MN.RecVal_Site)&!is.na(ModelParams_SIG$MN.RecVal_Sec))
  table(!is.na(ModelParams_SIG$MD.RecVal_Sec_All))
  
  #write.csv(RecSFMDataFrame,file = "data/Site_Recruits_SFM.csv")
  #write.csv(Jsec_P_Asec,file = "data/Sector_Recruits_REA.csv")
  #write.csv(ModelParams_SIG,file = "data/ModelParams_SIG.csv")
  
  # Two values of Site Level not matching!!!!
  
  #save(ModelParams_SIG, file = sprintf("data/%s_allmodfits_noNAs.rdata",Name))
  
  ModelParams_SIG %>% 
    group_by(Site,StartingYear) %>% 
    summarize(N=sum(g.N)) %>% 
    pivot_wider(names_from = StartingYear, values_from = N)
 
     
  
  ################################################################################
  ######### Regional + Genus models #################################
  ########## INTEGRATION ! #######################################################
  ################################################################################
  
  #Regional models (genus and region)
  combo_regional <- left_join(GrowthRegionalDF,RegionalAnnSurvDataFrame)
  combo_regional <- na.omit(combo_regional)
  
  #add SfM site-level stock recruitment data
  ModelParams_Regional = left_join(combo_regional,RecSFMSummary, by= c("Genus_Code","REGION") )
  
  #add REA sector-level stock recruitment data
  ModelParams_Regional = left_join(ModelParams_Regional,Regional_SectorStockRec,by=c("Genus_Code"="GENUS_CODE","REGION"))
  
  #save(ModelParams_Regional, file = sprintf("data/%s_RegionalModFits_noNAs.rdata",Name))
  

