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


# Loading Functions -------------------------------------------------------
leadz=function(x,n){return(formatC(x,width=n,flag=0))}
A2D=function(A){return(2*sqrt(A/pi))}
Am2Dcm=function(A){return(2*sqrt(A/pi)*100)}
D2A=function(D){return((D/2)^2*pi)}


# # Loading / Managing DataFrames: ColTrans, surv_dat  ----------------------------------------------  ---------
# ColTransAS=read.csv("./data/ColonyTransitions/ASRAMP23_ColonyTransitions.csv")
# ColTransMA=read.csv("./data/ColonyTransitions/MARAMP22_ColonyTransitions.csv")
# 
# ColTrans=rbind(ColTransAS,ColTransMA)
# ColTrans=ColTrans %>%
#   mutate(SiteInterval=paste0(ColTrans$Site,"_",ColTrans$Interval),
#          SArea_STAcm2=10^4*SArea_STA,
#          SArea_ENDcm2=10^4*SArea_END,
#          Shape_Area_STAcm2=10^4*Shape_Area_STA,
#          Shape_Area_ENDcm2=10^4*Shape_Area_END,
#          Shape_Leng_STAcm=10^2*Shape_Leng_STA,
#          Shape_Leng_ENDcm=10^2*Shape_Leng_END,
#   ) %>% 
#   rename(N_t0=nPatches_STA,
#          StartingDate=TL_Date_STA,
#          EndingDate=TL_Date_END,
#          TransitionType=TransitionTypeSimple
#   )
# 
# #Prepping surv_dat output
# survnames=c("Site_Genet","l10_Area_STA","Survival","Genus","Interval","SiteInterval","Site","N_t0",
#             "TransitionType","Interval_Years","StartingDate","EndingDate")
# surv_dat=ColTrans[,survnames]
# surv_dat=subset(surv_dat,TransitionType!="RECR")
# 
# # Change these file names when saving data to reflect the new changes
# # save(ColTrans, file="data/Colony_Data_20210917_edited.rdata")
# # save(surv_dat, file="data/Colony_Data_202100917_edited_survival.rdata")
# 
# 
# # Where we run the models -------------------------------------------------
# 
# #Define Site - Interval - Taxa groupings for which to run models
# 
# #SiteIntervalGenus
# ColTrans$SIG=paste0(substr(ColTrans$Site,5,8),substr(ColTrans$Site,9,11),"_",
#                     substr(year(ColTrans$StartingDate),3,4),leadz(month(ColTrans$StartingDate),2),"-",
#                     substr(year(ColTrans$EndingDate),3,4),leadz(month(ColTrans$EndingDate),2),"_",
#                     ColTrans$Genus)
# #SiteIntervalSpecies
# ColTrans$SIS=paste0(substr(ColTrans$Site,1,3),substr(ColTrans$Site,9,11),"_",
#                     substr(year(ColTrans$StartingDate),3,4),leadz(month(ColTrans$StartingDate),2),"-",
#                     substr(year(ColTrans$EndingDate),3,4),leadz(month(ColTrans$EndingDate),2),"_",
#                     ColTrans$Spec_Code)
# #Genus--> GENUS_CODE for simplicity later
# ColTrans=ColTrans %>% rename(GENUS_CODE=Genus)
# 
# #There are 7 entries without GENUS_CODE, dropping them for now
# length(which(is.na(ColTrans$GENUS_CODE)))
# GenetNAList=ColTrans[which(is.na(ColTrans$GENUS_CODE)),"Site_Genet"]
# ColTrans %>% filter(Site_Genet%in%GenetNAList) # All singleton genets (RECR/MORT)
# #Drop ʻEM - 11/18 I asked Corinne to figure these out...
# ColTrans=ColTrans %>% filter(!Site_Genet%in%GenetNAList)

# ColTransDup=read.csv("./Data/ColonyTransitions/FULLDOMAIN_19.22.23_ColonyTransitions.csv")
# dim(ColTransDup)
# dim(distinct(ColTransDup))
# ColTrans=distinct(ColTransDup)
# write.csv(ColTrans,"./Data/ColonyTransitions/FULLDOMAIN_19.22.23_ColonyTransitions_NoDup.csv")

ColTrans=read.csv("./Data/ColonyTransitions/FULLDOMAIN_19.22.23_ColonyTransitions_NoDup.csv")
ColTrans$REGION[which(ColTrans$REGION=="PRIA")]="PRIAs"
ColTrans$SEC_NAME[ColTrans$SEC_NAME%in%c("GUA_PITI_BOMB","GUA_TUMON_BAY")]="GUA_MP"
#Reassign OAH_XX_022 to OAH_OCC_005
ColTrans$Site[which(ColTrans$Site=="OAH_XXX_022")]="OAH_OCC_005"

#Report on Grouping
Usig=unique(ColTrans$SIG)
Usis=unique(ColTrans$SIS)
Nsig=length(Usig);Nsig
Nsis=length(Usis);Nsis

Name="HA_MA_AS_Models"

# Load Data to Get Proportion of "True" Recruitment ------------------------------------

#Assign SFM Sites to Sectors 
sitemd_MAAS=read.csv("./Data/MetaData/VitalRates_LatLong.csv");
sitemd_MAAS=sitemd_MAAS %>% rename(Site=ESD.Site.Name)

HAsitemd=read.csv("./Data/MetaData/HA_Site_Metadata.csv");
HAsitemd=HAsitemd %>% rename(Island=Island_Code,Site=SiteName) %>% 
  mutate(TP0=mdy(TP0),
         TP1=mdy(TP1),
         TP2=mdy(TP2),
         TP3=mdy(TP3),
         TP4=mdy(TP4),
         TP5=mdy(TP5),
         TP6=mdy(TP6)) %>% 
  pivot_longer(cols=TP0:TP6,values_to = "TimePointDate",names_to="TimePoint") %>% 
  na.omit() %>% 
  mutate(Year=year(TimePointDate)) %>% 
  select(Region,Island,Year,Site,Latitude,Longitude)
# names(sitemd)[1]="Site";names(sitemd)[1]="Site";sitemd$EndingDate=mdy(sitemd$SampleDate.MM.DD.YYYY.)
sitemd <- sitemd_MAAS %>% full_join(HAsitemd) %>% 
  arrange(Region,Island,Site) %>% 
  dplyr::select(Site, Latitude, Longitude) #get 1 lat/long for each site
sitemd <- distinct(sitemd)


uSD=ColTrans[,c("Site","StartingDate","EndingDate")] %>% 
  pivot_longer(cols=c("StartingDate","EndingDate"),values_to="Date") %>% 
  dplyr::select(all_of(c("Site","Date"))) %>% 
  unique()

CL.sf=as.data.frame(left_join(uSD,sitemd[,c("Site","Latitude","Longitude")],by=c("Site"))) #Site,Date,Lat,Long dataframe
CL.sf=st_as_sf(na.omit(CL.sf),coords=c("Longitude","Latitude")) #create geometry column
#secshp= st_read(dsn = "data/Shapefiles",layer = "ALLPacific_Sectors_Islands_5km_buffer")
secshp= st_read(dsn = "Data/Shapefiles",layer = "ALLPacific_Sectors_Islands")

#As: (1) we keep getting invalid geometries, and 
#    (2) local assignment is likely ok for planar geometry, we turn off spherical
# sf_use_s2(TRUE)
# secshp=st_make_valid(secshp)
# valid = st_is_valid(secshp)
# inval=which(!valid);print(inval)
sf_use_s2(FALSE)
secshp=st_make_valid(secshp)
valid = st_is_valid(secshp);inval=which(!valid);print(inval)

CL.sf = CL.sf %>% st_set_crs(value = st_crs(secshp)) #retrieve coordinate system
Site2Sec=st_join(CL.sf,secshp[,c("SEC_NAME","Region")])

#Add Sector to ColTrans dataframe!
Site2Sec_ONLY=unique(st_drop_geometry(Site2Sec)[,c("Site","SEC_NAME","Region")])
# Site2Sec_ONLY$SEC_NAME[Site2Sec_ONLY$Site=="FFS_OCC_002"]="French Frigate"
# Site2Sec_ONLY$SEC_NAME[Site2Sec_ONLY$Site=="FFS_OCC_014"]="French Frigate"
# Site2Sec_ONLY$SEC_NAME[Site2Sec_ONLY$Site=="HAW_OCC_003"]="HAW_PUNA"
# Site2Sec_ONLY$SEC_NAME[Site2Sec_ONLY$Site=="HAW_SIO_K08"]="HAW_KONA"

# ColTrans=ColTrans %>% left_join(Site2Sec_ONLY[,c("Site","SEC_NAME","Region")],by="Site")
# ColTrans=ColTrans %>% rename(REGION=Region)
# save(ColTrans, file="data/Colony_Data_MA_AS_20241112_edited.rdata")
Site2Sec_ONLY$SEC_NAME[Site2Sec_ONLY$SEC_NAME%in%c("GUA_PITI_BOMB","GUA_TUMON_BAY")]="GUA_MP"
Site2Sec_ONLY$SEC_NAME[Site2Sec_ONLY$Site%in%c("HAW_OCC_003")]="HAW_PUNA"
Site2Sec_ONLY$SEC_NAME[Site2Sec_ONLY$Site%in%c("HAW_SIO_K08")]="HAW_KONA"
#Get observed proportion of 'juveniles' as true recruits (SfM data)
#gray=juveniles, blue=true recruits
RECs=subset(ColTrans,TransitionType=="RECR")
bw=0.25
ggplot()+
  geom_histogram(data=subset(ColTrans,Shape_Area_ENDcm2>0),aes(A2D(A = Shape_Area_ENDcm2)),fill="grey",binwidth=bw)+
  xlim(c(-5*bw,100*bw))+
  geom_histogram(data=RECs,aes(A2D(Shape_Area_ENDcm2)),fill="blue",binwidth=bw)+
  facet_grid(GENUS_CODE~REGION,scales = "free_y")+
  geom_vline(xintercept = 5)+
  theme_bw()+
  ggtitle("Proportion of juveniles compared to true recruits")+
  labs(x="Diameter (cm)",color = "Legend")+scale_x_sqrt()
#Notes here: Area is in m2, but inter-survey interval is _large_ (i.e. 5 years, so "true" recruits are really adult cols)

#Check Proportion of Recruits per Taxon - hist calculation
#Model Mesh Points

binwidths=c(0.05,.1,.25,.5)#bw
uG=c("SSSS","ACSP","POCS","POSP","MOSP")
uR=unique(ColTrans$REGION)
PropRec_g=cbind(expand.grid(GENUS_CODE=uG,REGION=uR,BinWidth=binwidths),PropRec=NA)
up_breaks=5000
for(i_b in 1:length(binwidths)){
  binwidth=binwidths[i_b]
  N=5/binwidth
  for(i_R in 1:length(uR)){
    ColTransR=ColTrans %>% filter(REGION==uR[i_R],TransitionType%in%c("GROW","RECR"))
    RECsR=RECs%>% filter(REGION==uR[i_R])
    for(i_G in 1:length(uG)){
      if(uG[i_G]=="SSSS"){
        hcl=hist(A2D(subset(ColTransR,Shape_Area_ENDcm2>0)$Shape_Area_ENDcm2),breaks=seq(0,up_breaks,by=binwidth),plot = F)
        hrc=hist(A2D(RECsR$Shape_Area_ENDcm2),breaks=seq(0,up_breaks,by=binwidth),plot = F)
        PropRec_v=(hrc$counts[1:(N+1)]/hcl$counts[1:(N+1)])
        out_i=which(PropRec_g$GENUS_CODE==uG[i_G]&PropRec_g$REGION==uR[i_R]&PropRec_g$BinWidth==binwidth)
        PropRec_g$PropRec[out_i]=mean(PropRec_v[1:(N+1)],na.rm=T)
      }else{
        hcl=hist(A2D(subset(ColTransR,Shape_Area_ENDcm2>0&GENUS_CODE==uG[i_G])$Shape_Area_ENDcm2),breaks=seq(0,up_breaks,by=binwidth),plot = F)
        hrc=hist(A2D(subset(RECsR,GENUS_CODE==uG[i_G])$Shape_Area_ENDcm2),breaks=seq(0,up_breaks,by=binwidth),plot = F)
        PropRec_v=(hrc$counts[1:(N+1)]/hcl$counts[1:(N+1)])
        out_i=which(PropRec_g$GENUS_CODE==uG[i_G]&PropRec_g$REGION==uR[i_R]&PropRec_g$BinWidth==binwidth)
        PropRec_g$PropRec[out_i]=mean(PropRec_v[1:(N+1)],na.rm=T)
        # plot(hcl$breaks[1:(N+1)],PropRec_v,type="b",ylim=c(0,1),main=uG[i_G])
        # abline(h=mean(propREC,na.rm=T))}
      }
    }
  }
}
PropRec_g %>% arrange(GENUS_CODE,BinWidth)

#detach(package:plyr)
PropRecMn=PropRec_g %>% 
  group_by(GENUS_CODE,REGION) %>% 
  summarise(meanPropRec=mean(PropRec,na.rm = T))
#Final 0-5 cm diam. Proportion Recruits in the Juvenile size classes (will use for 'pro-rating')
as.data.frame(PropRecMn)



################################################################################
######### Site + Interval Years + Genus models #################################
########## Recruitment ! #######################################################
################################################################################

####Case 1 Assume Site - Level Stock Recruitment Relationship
####### Unlikely to hold true, but might be a decent approximation in skewed dispersal kernels
#recval is N recruits per Area of adult biomass - N/cm^2
# EQN 1: Then recval = ((Observed recruits-N)/ (Area Surveyed-cm^2)) / (Adult Area-cm^2)/Area Surveyed-cm^2)
# Adult area could come from site-level percent cover*area surveyed; or measured SFM adult area / Area surveyed

####Case 2 Assume Sector - Level Stock Recruitment Relationship

#recval is N recruits per Area of adult biomass - N/cm^2
# EQN 2: Then sector SFM recval = (SFM Observed recruits-N)/ Area Surveyed-cm^2) / (sector-level proportional cover) ###?
# EQN 3: Or sector REA recval = (REA Observed recruits***-N)/ Area Surveyed-cm^2) / (sector-level proportional cover)
# *** here we should prorate our REA juv density by SFM calibrated prop. of recruits in juv size classes


################################################################################
#SFM Observed Recruitment Rates 
################################################################################

#Case 1: Assume Site-Level Stock-Recruitment
# 
#Area Surveyed for each SIG
#Asurv=read.csv("./data/AreaSurveyed_N_Circrats.csv"); names(Asurv)[1]="Site"
#Asurv$Site[Asurv$Site=="OAH_XXX_022"]="OAH_OCC_005"
#Asurv$EndingDate=mdy(Asurv$Date)
# Asurv_l=Asurv[,c("Site","EndingDate","POCS","POSP","MOSP")]%>%
#   pivot_longer(cols=c("POSP","MOSP","POCS"),names_to=c("Genus"),values_to=c("Ncircrats")) %>%
#   mutate(Ncircrats=as.numeric(Ncircrats))
# Asurv_l$A.Surv.m2 =Asurv_l$Ncircrats*0.5 #numb of circrats * 0.5m2 (size of circrat) to get area surveyed in m2
# Asurv_l$A.Surv.cm2 = Asurv_l$Ncircrats*5000 #area surveyed cm2

Asurv_l_MAAS=read.csv("./data/MetaData/VitalRates_SurveyEffort.csv"); #names(Asurv)[1]="Site"
Asurv_l_MAAS=Asurv_l_MAAS %>% rename(A.Surv.m2=Effort,GENUS_CODE=Genus) %>% mutate(A.Surv.cm2=10^4*A.Surv.m2)

Asurv_lHA=read.csv("./data/MetaData/HA_AreaSurveyed_N_Circrats.csv"); #names(Asurv)[1]="Site"
Asurv_lHA$Site[Asurv_lHA$Site=="OAH_XXX_022"]="OAH_OCC_005"
Asurv_lHA=Asurv_lHA %>%
  pivot_longer(cols = c("POCS","POSP","MOSP"),names_to = "GENUS_CODE",values_to = "N.quads") %>% 
  mutate(Island=substr(Site,1,3),
         N.quads=as.numeric(N.quads),
         A.Surv.m2=N.quads*0.5,
         A.Surv.cm2=10^4*A.Surv.m2) %>% 
  select(c(Island,Site,GENUS_CODE,A.Surv.m2,A.Surv.cm2))
Asurv_l=full_join(Asurv_l_MAAS,Asurv_lHA) %>% arrange(Site) %>% distinct()
#table(Asurv_l$Site,Asurv_l$GENUS_CODE)
Asurv_l=Asurv_l %>%
  group_by(Island,Site,GENUS_CODE) %>% 
  summarize(A.Surv.m2=max(A.Surv.m2),
            A.Surv.cm2=max(A.Surv.cm2))
#table(Asurv_l$Site,Asurv_l$GENUS_CODE)

#Area of Adult Colonies for each SIG
Atax=ColTrans %>%
  group_by(SIG,Site,Interval,GENUS_CODE,StartingDate,EndingDate,Interval_Years) %>%
  summarise(AdultCoralArea_cm2=sum(Shape_Area_STAcm2)) #Area of adult cm2 = sum all corals in each genus for each year

#Link to Nrecruits from the Shape_Area_STAcm2 summed area: Adults at beginning of interval generate
RecSFMTib=ColTrans %>%
  filter(TransitionType %in% c("RECR")) %>%
  group_by(REGION,Island,SEC_NAME,Site,Interval,GENUS_CODE,SIG,StartingDate,EndingDate,Interval_Years) %>%
  summarize(Nrec=length(which(TransitionType=="RECR"))) %>%
  left_join(Atax,by=c("Site","Interval","GENUS_CODE","SIG","StartingDate","EndingDate","Interval_Years")) %>%
  left_join(Asurv_l,by=c("Island","Site","GENUS_CODE")) %>%
  mutate(
    # numb recruits/area surveyed cm2
    RecSFM_p_Survcm2=Nrec/A.Surv.cm2,
    #Case 1: Assume Site-Level Stock-Recruitment
    #Recruit N per adult area and annual rate
    RecSFM_p_SiteAdcm2=Nrec/AdultCoralArea_cm2,
    RecSFM_p_SiteAdcm2_Yr=RecSFM_p_SiteAdcm2/Interval_Years) # numb recruits/total adult area cm2 annualized

RecSFMDataFrame=as.data.frame(RecSFMTib)

save(RecSFMDataFrame, file = paste0("./Data/ModelData/",Name,"_RecSFMrates.rdata"))

#plot
hist(RecSFMDataFrame$RecSFM_p_SiteAdcm2_Yr, breaks = 20)


#REGIONAL CALCULATIONS BY GENUS AND REGION
RecSFMDataFrame$RecSFM_p_SiteAdcm2_Yr[RecSFMDataFrame$RecSFM_p_SiteAdcm2_Yr == Inf] <-NA
#Calculate median, mean, and lower/upper quantiles for site-level stock recruitment GENUS values
RecSFMSummary=RecSFMDataFrame %>%
  dplyr::select(SEC_NAME,Site,Interval,GENUS_CODE,REGION,RecSFM_p_SiteAdcm2_Yr)%>%
  group_by(GENUS_CODE,REGION)%>%
  summarise(SiteStock_median = median(RecSFM_p_SiteAdcm2_Yr, na.rm = TRUE),
            SiteStock_mean = mean(RecSFM_p_SiteAdcm2_Yr, na.rm = TRUE),
            SiteStock_q05= quantile(RecSFM_p_SiteAdcm2_Yr, c(0.05),na.rm = TRUE),
            SiteStock_q95= quantile(RecSFM_p_SiteAdcm2_Yr, c(0.95),na.rm = TRUE)
  )



#REA Observed Recruitment Rates per adult cover, at Sector and Site Level,
#with SIG matching site-level data grouped at sector containing SIG

#Case 2: Assum Sector-Level Stock-Recruitment
#Site Recruitment, Sector Stock
#Get Site-Level and Sector_level REA juv density
#rea_sec=read.csv("data/NWHI_MHI_REA_Data_Sector.csv")

#Get Percent Cover by Genus at each Sector
# cov1_sec=read.csv("../data/NWHI_MHI_Cover_T1_Data_Sector.csv")
# cov3_sec=read.csv("data/NWHI_MHI_Cover_T3_Data_Sector.csv")
uSec=na.omit(unique(ColTrans$SEC_NAME))
cov1_sec=read.csv("T:/Benthic/Data/REA Coral Demography & Cover/Summary Data/Sector/BenthicCover_2010-2023_Tier1_SECTOR_updated.csv")
names(cov1_sec)[1:16]=substr(names(cov1_sec)[1:16],6,999)
meta=names(cov1_sec)[1:5]
classes1=names(cov1_sec)[6:14]
SEclasses1=names(cov1_sec)[22:30]

cov3_sec=read.csv("T:/Benthic/Data/REA Coral Demography & Cover/Summary Data/Sector/BenthicCover_2010-2023_Tier3_SECTOR_updated.csv")
names(cov3_sec)[1:110]=substr(names(cov3_sec)[1:110],6,999)
classes3=names(cov3_sec)[6:110]
SEclasses3=names(cov3_sec)[117:(length(classes3)+117-1)]

#select first 15 cols of cov1_sec and pivot CORAL,MA,TURF,CCA,EMA,SC,I,SED,HAL to
#a new column titled CAT and values to column titled COVER
cov1_lse = cov1_sec[,c(meta,SEclasses1)] %>% 
  pivot_longer(cols = all_of(SEclasses1),names_to="CAT",values_to="SE.COVER")
cov1_lse$CAT=substr(cov1_lse$CAT,10,999) #drop pooled SE
cov1_l = cov1_sec[,c(meta,classes1)] %>% #get cover by category & N
  pivot_longer(cols = all_of(classes1),names_to="CAT",values_to="COVER") %>% 
  left_join(cov1_lse,by=c(meta,"CAT"))
cov1_l=cov1_l[,c("REGION","ISLAND","ANALYSIS_SEC","ANALYSIS_YEAR","CAT","COVER","SE.COVER","N")]
cov1_l$SD.COVER=cov1_l$SE.COVER*sqrt(cov1_l$N)
cov1_l$VAR.COVER=cov1_l$SD.COVER^2

cov3_lse = cov3_sec[,c(meta,SEclasses3)] %>% 
  pivot_longer(cols = all_of(SEclasses3),names_to="CAT",values_to="SE.COVER")
cov3_lse$CAT=substr(cov3_lse$CAT,10,999)
cov3_l = cov3_sec[,c(meta,classes3)] %>% 
  pivot_longer(cols = all_of(classes3),names_to="CAT",values_to="COVER") %>% 
  left_join(cov3_lse,by=c(meta,"CAT"))
cov3_l=cov3_l[,c("REGION","ISLAND","ANALYSIS_SEC","ANALYSIS_YEAR","CAT","COVER","SE.COVER","N")]
cov3_l$SD.COVER=cov3_l$SE.COVER*sqrt(cov3_l$N)
cov3_l$VAR.COVER=cov3_l$SD.COVER^2


#including all acropora morphs to get best recruit #s
MPPcovers=c("ACBR","ACTA","MOBR","MOEN","MONE","MOFO","POCS","POEN","POFO","POMA","PONM","POBR") #cover categories for each genus for Tier3 benthic cat
names(MPPcovers)=c(rep("ACSP",2),rep("MOSP",4),"POCS",rep("POSP",5))
MMPlu=names(MPPcovers)
names(MMPlu) = MPPcovers

cov_l=rbind(subset(cov1_l,CAT=="CORAL"),subset(cov3_l,CAT%in%MPPcovers)) 
cov_l$GENUS_CODE=MMPlu[cov_l$CAT] #add Genus column
cov_l$GENUS_CODE[cov_l$CAT=="CORAL"]="SSSS" #if category is coral, assign SSSS genus code
cov_l$SEC_YEAR=paste0(cov_l$ANALYSIS_SEC,"-",cov_l$ANALYSIS_YEAR)
cov_l$REGION=factor(cov_l$REGION,levels=c("MHI","NWHI","PRIAs","MARIAN","SAMOA"))
#SUM all within genus variation at each site
keepcols=c("REGION","ISLAND","ANALYSIS_SEC","ANALYSIS_YEAR","GENUS_CODE","SEC_YEAR")
cov_sum=cov_l %>%  #cover by genus code and sector
  group_by(ANALYSIS_SEC,ANALYSIS_YEAR,GENUS_CODE) %>% 
  summarize(COVER=sum(COVER),VAR.COVER=sum(VAR.COVER),N.COVER=sum(N))
cov_sum$SD.COVER=sqrt(cov_sum$VAR.COVER)
cov_sum$SE.COVER=cov_sum$SD.COVER/sqrt(cov_sum$N.COVER)

cov=left_join(cov_sum,unique(cov_l[,keepcols]),by=c("ANALYSIS_SEC","ANALYSIS_YEAR","GENUS_CODE")) #cov_l by distint benth cats. sum all cover by genus
cov=cov[,c(keepcols,"COVER","SE.COVER","N.COVER")] #reorganize columns
#cover is percent cover (0-100%)


# prorated mean juvenile col density for each sector and each year divided by adult cover sector year genus
#get mean, 5 and 95 quantiles for each sector year genus
#for each of those ^^ 3 numbers, want to prorate depending on genus code. 
ISL_luDF=cov %>% ungroup() %>% select(ANALYSIS_SEC,ISLAND) %>% distinct()
rea_sec=read.csv("T:/Benthic/Data/REA Coral Demography & Cover/Summary Data/Sector/BenthicREA_sectordata_GENUS.csv")#read.csv("../fish-paste/data/Sectors-Strata-Areas.csv")
rea_sec=rea_sec %>%
  rename(ANALYSIS_SEC=PooledSector_Demo_Viztool,
         N.JCD=n)
#setdiff(ISL_luDF$ANALYSIS_SEC,rea_sec$ANALYSIS_SEC)
ISL_luDF=rbind(ISL_luDF,data.frame(ANALYSIS_SEC=setdiff(rea_sec$ANALYSIS_SEC,ISL_luDF$ANALYSIS_SEC),ISLAND=c("Swains","Tau","Tutuila","Tutuila")))
rea_sec=rea_sec %>%
  left_join(ISL_luDF)
head(rea_sec)
#which(is.na(rea_sec$ISLAND))

##############################################
#NEEDS SEC_NAME REVERSION for MHI
##############################################
table(cov$ANALYSIS_YEAR,useNA = "always")
table(rea_sec$ANALYSIS_YEAR,useNA = "always")
names(rea_sec)[which(names(rea_sec)=="n")]="N.JCD"
#rea_sec$OBS_YEAR=rea_sec$ANALYSIS_YEAR

# Propagating error for sector level stock recruitment
#mJCDp =  mJCD*MeanPropRec
#seJCDp = seJCD*MeanPropRec
#sdJCDp = sqrt(nJCD)*seJCDp
#nJCDp = nJCD
#sdAC = sqrt(nAC)*seAC
#mRec = mJCDp/mAC
#sdRec = mRec*sqrt((sdJCDp/mJCDp)^2 + (sdAC/mAC)^2))
#nRec = min(nJCD,nAC)????
#seRec = sdRec/sqrt(nRec)
#CI95Rec = 1.96*seRec

cov$P_COVER=cov$COVER/100 #convert percent cover to proportional cover
cov$SE.P_COVER=cov$SE.COVER/100
rea_sec$Mean_JuvColDen_cm2=  rea_sec$Mean_JuvColDen/10^4 #convert mean juv col density to cm2
rea_sec$SE_JuvColDen_cm2=  rea_sec$SE_JuvColDen/10^4 #convert juv col den SE to cm2
#Sector Mismatch Issues Between ColTrans for TAU_OPEN and TUT_FAGATELE
# ColTrans %>% filter(SEC_NAME %in%c("TAU_OPEN","TUT_FAGATELE"))%>% nrow()
# ColTrans %>% filter(SEC_NAME%in%c("TAU_ALL","TUT_FAGATELE_FAGALUA"))%>% nrow()
# cov %>% filter(ANALYSIS_SEC%in%c("TAU_OPEN","TUT_FAGATELE")) %>% nrow()
# cov %>% filter(ANALYSIS_SEC%in%c("TAU_ALL","TUT_FAGATELE_FAGALUA"))%>% nrow()
# rea_sec %>% filter(ANALYSIS_SEC%in%c("TAU_OPEN","TUT_FAGATELE"))%>% nrow()
# rea_sec %>% filter(ANALYSIS_SEC%in%c("TAU_ALL","TUT_FAGATELE_FAGALUA"))%>% nrow()
rea_sec$ANALYSIS_SEC[which(rea_sec$ANALYSIS_SEC=="TAU_ALL")]="TAU_OPEN"
rea_sec$ANALYSIS_SEC[which(rea_sec$ANALYSIS_SEC=="TUT_FAGALUA_FAGATELE")]="TUT_FAGATELE"

#Also REGION Mismatches
sort(unique(ColTrans$REGION))
sort(unique(as.character(cov$REGION)))
sort(unique(rea_sec$REGION))
rea_sec$REGION[which(rea_sec$REGION=="CNMI")]="MARIAN"
rea_sec$REGION[which(rea_sec$REGION=="GUA")]="MARIAN"

#Also need to better match up "ANALYSIS_YEAR"S between cover and rea
rea_sec$ANALYSIS_YEAR[rea_sec$ANALYSIS_YEAR=="2013"]="2013-15"
#Calculate sector level stock recruitment and propagate error
#sector level recruitment for sectors for this study
cov$ANALYSIS_YEAR[which(cov$ISLAND%in%c("Baker","Howland")&cov$ANALYSIS_YEAR=="2015")]="2014-15"
cov$ANALYSIS_YEAR[which(cov$ISLAND%in%c("Baker","Howland")&cov$ANALYSIS_YEAR=="2017")]="2017-18"

#
cov$ANALYSIS_SEC[which(cov$ANALYSIS_SEC%in%c("TUT_NE"))]="TUT_NE_OPEN"
cov$ANALYSIS_SEC[which(cov$ANALYSIS_SEC%in%c("TUT_NW"))]="TUT_NW_OPEN"
cov$ANALYSIS_SEC[which(cov$ANALYSIS_SEC%in%c("TUT_SE"))]="TUT_SE_OPEN"
cov$ANALYSIS_SEC[which(cov$ANALYSIS_SEC%in%c("TUT_AUNUU_B"))]="TUT_AUNUU"
rea_sec$ANALYSIS_YEAR[which(rea_sec$REGION%in%c("SAMOA")&rea_sec$ANALYSIS_YEAR=="2015")]="2015-16"
cov$ANALYSIS_SEC[which(cov$ANALYSIS_SEC=="GUA_MP_MINUS_ACHANG")]="GUA_MP"
########################################################## Stopped Here 11/18 1100
names(cov)[which(names(cov)=="Genus")]="GENUS_CODE"
Jsec_P_Asec= rea_sec %>%
  filter(ANALYSIS_SEC%in%uSec&GENUS_CODE%in%c("SSSS","ACSP","MOSP","POSP","POCS")) %>% #filter by target taxa and target sectors (uSec is unique sector list from above)
  select(REGION,ISLAND,ANALYSIS_SEC,ANALYSIS_YEAR,GENUS_CODE,Mean_JuvColDen_cm2, SE_JuvColDen_cm2,N.JCD)%>%
  left_join(cov,by=c("REGION","ISLAND","ANALYSIS_SEC","GENUS_CODE","ANALYSIS_YEAR")) %>% 
  left_join(PropRecMn,by=c("REGION","GENUS_CODE")) %>%
  group_by(REGION,ISLAND,ANALYSIS_SEC,ANALYSIS_YEAR,GENUS_CODE)%>%
  mutate(Mean_JuvColDenP=Mean_JuvColDen_cm2*meanPropRec,
         SE_JuvColDenP=SE_JuvColDen_cm2*meanPropRec,
         N.JCDp=N.JCD,
         SD_JuvColDenP=sqrt(N.JCDp)*SE_JuvColDenP,
         SD.P_COVER=sqrt(N.COVER)*SE.P_COVER,
         MN.RecVal = Mean_JuvColDenP/P_COVER,
         SE.RecVal = MN.RecVal*sqrt((SE_JuvColDenP/Mean_JuvColDenP)^2+(SE.P_COVER/P_COVER)^2),#propogating SE error
         N.RecVal = min(N.JCDp,N.COVER,na.rm=T),
         #SE.RecVal = SD.RecVal/sqrt(N.RecVal),
         CI95.RecVal = 1.96*SE.RecVal,
         LOW_CI95.RecVal= MN.RecVal - CI95.RecVal,
         HIGH_CI95.RecVal= MN.RecVal + CI95.RecVal
  )

#sector level recruitment for all sectors in the NWHI and MHI (not just study sectors)
Jsec_P_ALL= rea_sec %>%
  filter(GENUS_CODE%in%c("SSSS","ACSP","MOSP","POSP","POCS")) %>% #uSec is unique sector list from above
  select(REGION,ISLAND,ANALYSIS_SEC,ANALYSIS_YEAR,GENUS_CODE,Mean_JuvColDen_cm2, SE_JuvColDen_cm2,N.JCD)%>%
  left_join(cov,by=c("REGION","ISLAND","ANALYSIS_SEC","GENUS_CODE","ANALYSIS_YEAR")) %>% 
  left_join(PropRecMn,by=c("REGION","GENUS_CODE")) %>%
  group_by(REGION,ISLAND,ANALYSIS_SEC,ANALYSIS_YEAR,GENUS_CODE)%>%
  mutate(Mean_JuvColDenP=Mean_JuvColDen_cm2*meanPropRec,
         SE_JuvColDenP=SE_JuvColDen_cm2*meanPropRec,
         N.JCDp=N.JCD,
         SD_JuvColDenP=sqrt(N.JCDp)*SE_JuvColDenP,
         SD.P_COVER=sqrt(N.COVER)*SE.P_COVER,
         MN.RecVal = Mean_JuvColDenP/P_COVER,
         SE.RecVal = MN.RecVal*sqrt((SE_JuvColDenP/Mean_JuvColDenP)^2+(SE.P_COVER/P_COVER)^2),#propogating SE error
         N.RecVal = min(N.JCDp,N.COVER,na.rm=T),
         #SE.RecVal = SD.RecVal/sqrt(N.RecVal),
         CI95.RecVal = 1.96*SE.RecVal,
         LOW_CI95.RecVal= MN.RecVal - CI95.RecVal,
         HIGH_CI95.RecVal= MN.RecVal + CI95.RecVal
  )

Jsec_P_ALL %>% filter(ISLAND%in%Jsec_P_Asec$ISLAND) %>% #c("Hawaii","Maui","Oahu","French Frigate","Lisianski","Kure")
  ggplot(aes(MN.RecVal))+geom_histogram(binwidth=0.01)+facet_grid(ISLAND~GENUS_CODE,scales="free")
#View(Jsec_P_Asec)  

#rec values by region and genus code
RecVal_Sec_Dists=Jsec_P_ALL %>% filter(!is.infinite(MN.RecVal)) %>% 
  group_by(REGION,GENUS_CODE) %>% 
  summarize(MD.RecVal_Sec_All=median(MN.RecVal,na.rm=T),
            MD.CI95_LO.RecVal_Sec_All=median(MN.RecVal-1.96*SE.RecVal,na.rm=T),
            MD.CI95_HI.RecVal_Sec_All=median(MN.RecVal+1.96*SE.RecVal,na.rm=T),
            MN.RecVal_Sec_All=mean(MN.RecVal,na.rm=T),
            MN.CI95_LO.RecVal_Sec_All=mean(MN.RecVal-1.96*SE.RecVal,na.rm=T),
            MN.CI95_HI.RecVal_Sec_All=mean(MN.RecVal+1.96*SE.RecVal,na.rm=T)
  )

Jsec_P_ALL %>% filter(ISLAND%in%Jsec_P_Asec$ISLAND) %>% #filter(ISLAND%in%c("Hawaii","Maui","Oahu","French Frigate","Lisianski","Kure")) %>% 
  ggplot(aes(MN.RecVal))+
  geom_histogram(binwidth=0.01)+
  geom_vline(aes(xintercept=MD.RecVal_Sec_All),data=RecVal_Sec_Dists,color="darkcyan",lty=1)+
  geom_vline(aes(xintercept=MD.CI95_LO.RecVal_Sec_All),data=RecVal_Sec_Dists,color="darkcyan",lty=3)+
  geom_vline(aes(xintercept=MD.CI95_HI.RecVal_Sec_All),data=RecVal_Sec_Dists,color="darkcyan",lty=3)+
  facet_grid(REGION~GENUS_CODE)+
  scale_x_sqrt()

#plot site-level stock recruitment values and sector-level stock recruitment values
bw=0.01
recvals=c(Jsec_P_ALL$MN.RecVal,RecSFMDataFrame$RecSFM_p_SiteAdcm2_Yr);recvals[is.infinite(recvals)]=NA
xbnds=c(-2*bw,max(recvals,na.rm=T))
Site=ggplot(RecSFMDataFrame,aes(RecSFM_p_SiteAdcm2_Yr))+
  geom_histogram(binwidth=bw)+ theme_bw()+
  facet_grid("GENUS_CODE")+
  #facet_wrap(Genus~REGION)+
  xlim(xbnds)+scale_x_sqrt()+
  xlab("Site-level stock recruitment") # "# recruits  / area adult coral"

Sec= Jsec_P_Asec%>%
  filter(GENUS_CODE%in%c("ACSP","MOSP","POSP","POCS")) %>%
  ggplot(aes(MN.RecVal))+
  geom_histogram(binwidth=bw)+theme_bw()+
  facet_grid("GENUS_CODE")+
  #facet_wrap(Genus~REGION)+
  xlim(xbnds)+scale_x_sqrt()+
  xlab("Proportional sector-level stock recruitment") # "Proportional juv. colony density / cover"

Site/Sec

#Sector Stock Rec for REGIONAL MODEL
Regional_SectorStockRec <- as.data.frame(Jsec_P_Asec) 
Regional_SectorStockRec$MN.RecVal[Regional_SectorStockRec$MN.RecVal == Inf] <- NA
Regional_SectorStockRec$SE.RecVal[Regional_SectorStockRec$SE.RecVal == Inf] <- NA

Regional_SectorStockRec <- Regional_SectorStockRec %>%
  select(ANALYSIS_SEC,ANALYSIS_YEAR,GENUS_CODE,REGION,MN.RecVal,SE.RecVal)%>%
  group_by(GENUS_CODE, REGION) %>%
  summarise(SectorStock_median = median(MN.RecVal,na.rm = TRUE),
            SectorStock_MD_CI95_LO = median(MN.RecVal-1.96*SE.RecVal,na.rm=T),
            SectorStock_MD_CI95_HI = median(MN.RecVal+1.96*SE.RecVal,na.rm=T),
            SectorStock_mean = mean(MN.RecVal,na.rm = TRUE),
            SectorStock_MN_CI95_LO = mean(MN.RecVal-1.96*SE.RecVal, na.rm = TRUE), 
            SectorStock_MN_CI95_HI = mean(MN.RecVal+1.96*SE.RecVal, na.rm = TRUE) 
  )
#change negative low CI values to 0 since you can't have negative recruitment
Regional_SectorStockRec[Regional_SectorStockRec<0] = 0

#Include Sec_All values in Regional sector stock
Regional_SectorStockRec <- left_join(Regional_SectorStockRec,RecVal_Sec_Dists, by= c("GENUS_CODE","REGION"))

#Check unit consistency of log size attributes
sort(unique(ColTrans$TransitionType))
ColTrans$TransitionType[ColTrans$TransitionType=="SHRINK"]="GROW"
inf2NA=function(x){x[which(is.infinite(x))]=NA;return(x)}
ColTrans=ColTrans %>% mutate(
  log10_SS=inf2NA(log10(Shape_Area_STAcm2)),
  log10_ES=inf2NA(log10(Shape_Area_ENDcm2)),               
  l10TransitionMagnitude=log10_ES-log10_SS,
  log10_ESc=log10_SS+(l10TransitionMagnitude/Interval_Years),
  log10_SSsa=inf2NA(log10(SArea_STAcm2)),
  log10_ESsa=inf2NA(log10(SArea_ENDcm2)),
  l10STransitionMagnitude=log10_ESsa-log10_SSsa,
  log10_EScsa=log10_SSsa+(l10STransitionMagnitude/Interval_Years),
  log10(Shape_Area_STAcm2),Circularity=(4*pi*Shape_Area_STA)/(Shape_Leng_STA^2),
  NfragBin=cut(N_t0,breaks=c(0,1,3,10,100)+.5,labels=c("1","2-3","4-9","10+")),
  CircularityBin=cut(Circularity,breaks=seq(0,1,length.out = 5)),
  AnnualPropRate_E=1+(10^(log10_ESc)-10^(log10_SS))/10^log10_SS,
  DailyPropRate_E=AnnualPropRate_E^(1/365),
  MonthlyPropRate_E=AnnualPropRate_E^(1/12)
)

#Output Data.Frames to continue
save(list=c("Site2Sec","Site2Sec_ONLY","ColTrans","RecSFMDataFrame","RecSFMSummary","Regional_SectorStockRec","Jsec_P_Asec","RecVal_Sec_Dists","Jsec_P_ALL"),
     file = "./Data/ColonyTransitions/Script_Step3_DataPackage.rdata")


#Sector check
all(sort(unique(ColTrans$SEC_NAME))==sort(unique(Site2Sec_ONLY$SEC_NAME)))
all(sort(unique(ColTrans$SEC_NAME))==sort(unique(RecSFMDataFrame$SEC_NAME)))
all(sort(unique(ColTrans$SEC_NAME))==sort(unique(Jsec_P_Asec$ANALYSIS_SEC))) # Don't have juvenile counts for TUT_FAGATELE and TAU_OPEN
View(ColTrans)
