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
D2A=function(D){return((D/2)^2*pi)}
source("R/gcdist.R")

# Loading / Managing DataFrames: ColonyLevel, surv_dat  ----------------------------------------------  ---------
# Already Done?
load("data/Colony_Data_20210917_edited.rdata") #ColonyLevel
#load("data/Colony_Data_202100917_edited_survival.rdata") #surv_dat
#skip to line 60!

# If not, edit ColonyLevel here!
##Manage ColonyLevel and create surv_dat
load("data/Patch_And_Colony_Data_20210917.rdata") #ColonyLevel
if(!all(c("Survival","log10_ESc")%in%names(ColonyLevel))){
  #LN size
  ColonyLevel$ln_SS <- log(ColonyLevel$StartingSize)
  ColonyLevel$ln_ES <- log(ColonyLevel$EndingSize)
  #Log10 size
  ColonyLevel$log10_SS <- log10(ColonyLevel$StartingSize)
  ColonyLevel$log10_ES <- log10(ColonyLevel$EndingSize)
  #changing ln_ES/ln_SS from Na/-Inf to an actual value for recruitment/mortality
  ColonyLevel$log10_SS <- as.numeric(ColonyLevel$log10_SS)
  ColonyLevel$log10_SS[ColonyLevel$log10_SS == -Inf] <-NA
  
  #Mortality to Survival
  ColonyLevel <- cbind(ColonyLevel, data.frame(Survival = 1 - ColonyLevel$Mortality))
  
  # Making the size term yearly
  GY <- (ColonyLevel$ln_ES - ColonyLevel$ln_SS)*(1/ColonyLevel$Interval_Year)
  ColonyLevel <- cbind(ColonyLevel, data.frame(ln_ESc = ColonyLevel$ln_SS + GY))
  GY <- (ColonyLevel$log10_ES - ColonyLevel$log10_SS)*(1/ColonyLevel$Interval_Year)
  ColonyLevel <- cbind(ColonyLevel, data.frame(log10_ESc = ColonyLevel$log10_SS + GY))
  ColonyLevel$log10_ESc[ColonyLevel$log10_ESc == -Inf] <-NA
}
#Add island code to Lisianski
ColonyLevel$Island <- paste0(substr(ColonyLevel$Site,1,3))
#add region to ColonyLevel datset
RegionLU=c("NWHI","MHI","NWHI","NWHI","MHI","MHI","NWHI") #create a lookup table
names(RegionLU)=c("FFS","HAW","KUR","LIS","MAI","OAH","PHR")
ColonyLevel$REGION = RegionLU[ColonyLevel$Island]


#Prepping surv_dat output
surv_dat=ColonyLevel[,c("ColonyID","log10_SS","Survival","Genus_Code","Interval","SiteInterval","Site","N_t0","TransitionType","Interval_Years","StartingDate","EndingDate")]
surv_dat=subset(surv_dat,TransitionType!="RECR")
names(surv_dat)=c("ColonyID","size","survival","Genus_Code","Interval","SiteInterval","Site","N_t0","TransitionType","Interval_Years")

# Change these file names when saving data to refect the new changes
# save(ColonyLevel, file="data/Colony_Data_20210917_edited.rdata")
# save(surv_dat, file="data/Colony_Data_202100917_edited_survival.rdata")


# Where we run the models -------------------------------------------------

#Define Site - Interval - Taxa groupings for which to run models
#SiteIntervalGenus
ColonyLevel$SIG=paste0(substr(ColonyLevel$Site,1,3),substr(ColonyLevel$Site,9,11),"_",
                       substr(year(ColonyLevel$StartingDate),3,4),leadz(month(ColonyLevel$StartingDate),2),"-",
                       substr(year(ColonyLevel$EndingDate),3,4),leadz(month(ColonyLevel$EndingDate),2),"_",
                       ColonyLevel$Genus_Code)
#SiteIntervalSpecies
ColonyLevel$SIS=paste0(substr(ColonyLevel$Site,1,3),substr(ColonyLevel$Site,9,11),"_",
                       substr(year(ColonyLevel$StartingDate),3,4),leadz(month(ColonyLevel$StartingDate),2),"-",
                       substr(year(ColonyLevel$EndingDate),3,4),leadz(month(ColonyLevel$EndingDate),2),"_",
                       ColonyLevel$Spec_Code)

#Report on Grouping
Usig=unique(ColonyLevel$SIG)
Usis=unique(ColonyLevel$SIS)
Nsig=length(Usig);Nsig
Nsis=length(Usis);Nsis

Name="Streamlined_VR_Models_"

# Load Data to Get Proportion of "True" Recruitment ------------------------------------
#Assign SFM Sites to Sectors 
sitemd=read.csv("data/Thesis_Sites_Metadata2.csv");names(sitemd)[1]="Site";names(sitemd)[1]="Site";sitemd$EndingDate=mdy(sitemd$SampleDate.MM.DD.YYYY.)
sitemd <- sitemd %>% dplyr::select(Site, Latitude, Longitude) #get 1 lat/long for each site
sitemd <- distinct(sitemd)
#sitemd=sitemd %>% group_by(Site) %>% summarize(Lat=mean(Latitude),Lon=mean(Longitude))

#Reassign OAH_XX_022 to OAH_OCC_005
ColonyLevel$Site[which(ColonyLevel$Site=="OAH_XXX_022")]="OAH_OCC_005"

uSD=ColonyLevel[,c("Site","StartingDate","EndingDate")] %>% 
  pivot_longer(cols=c("StartingDate","EndingDate"),values_to="Date") %>% 
  dplyr::select(all_of(c("Site","Date"))) %>% 
  unique()

CL.sf=as.data.frame(left_join(uSD,sitemd[,c("Site","Latitude","Longitude")],by=c("Site"))) #Site,Date,Lat,Long dataframe
CL.sf=st_as_sf(na.omit(CL.sf),coords=c("Longitude","Latitude")) #create geometry column
secshp= st_read(dsn = "data/SectorSHP",layer = "ALLPacific_Sectors_Islands_5km_buffer")
if(!all(st_is_valid(secshp))){secshp=st_make_valid(secshp)}
st_is_valid(secshp)[-26]

CL.sf = CL.sf %>% st_set_crs(value = st_crs(secshp[-26])) #retrieve coordinate system
Site2Sec=st_join(CL.sf,secshp[-26,"SEC_NAME"])

#Add Sector to ColonyLevel dataframe!
Site2Sec_ONLY=unique(st_drop_geometry(Site2Sec)[,c("Site","SEC_NAME")])
Site2Sec_ONLY$SEC_NAME[Site2Sec_ONLY$Site=="FFS_OCC_002"]="French Frigate"
Site2Sec_ONLY$SEC_NAME[Site2Sec_ONLY$Site=="FFS_OCC_014"]="French Frigate"
Site2Sec_ONLY$SEC_NAME[Site2Sec_ONLY$Site=="HAW_OCC_003"]="HAW_PUNA"
Site2Sec_ONLY$SEC_NAME[Site2Sec_ONLY$Site=="HAW_SIO_K08"]="HAW_KONA"

#ColonyLevel=ColonyLevel %>% left_join(Site2Sec_ONLY[,c("Site","SEC_NAME")],by="Site")
#save(ColonyLevel, file="data/Colony_Data_20210917_edited.rdata")


#Get observed proportion of 'juveniles' as true recruits (SfM data)
#gray=juveniles, blue=true recruits
RECs=subset(ColonyLevel,TransitionTypeSimple=="RECR")
bw=0.5
ggplot()+
  geom_histogram(data=subset(ColonyLevel,EndingSize>0),aes(A2D(EndingSize)),fill="grey",binwidth=bw)+
  xlim(c(0,10))+
  geom_histogram(data=RECs,aes(A2D(EndingSize)),fill="blue",binwidth=bw)+
  facet_grid(Genus_Code~.,scales = "free_y")+
  geom_vline(xintercept = 5)+
  theme_bw()+
  ggtitle("Proportion of juveniles compared to true recruits")+
  labs(x="Diameter (cm)",color = "Legend")

#Check Proportion of Recruits per Taxon - hist calculation
binwidths=bw
uG=c("SSSS","POCS","POSP","MOSP")
PropRec_g=cbind(expand.grid(Genus_Code=uG,BinWidth=binwidths),PropRec=NA)
for(i_b in 1:length(binwidths)){
  binwidth=binwidths[i_b]
  N=5/binwidth
  for(i_G in 1:length(uG)){
    if(uG[i_G]=="SSSS"){
      hcl=hist(A2D(subset(ColonyLevel,EndingSize>0)$EndingSize),breaks=seq(0,1000,by=binwidth),plot = F)
      hrc=hist(A2D(RECs$EndingSize),breaks=seq(0,1000,by=binwidth),plot = F)
      PropRec_v=(hrc$counts[1:(N+1)]/hcl$counts[1:(N+1)])
      out_i=which(PropRec_g$Genus_Code==uG[i_G]&PropRec_g$BinWidth==binwidth)
      PropRec_g$PropRec[out_i]=mean(PropRec_v[1:(N+1)],na.rm=T)
    }else{
      hcl=hist(A2D(subset(ColonyLevel,EndingSize>0&Genus_Code==uG[i_G])$EndingSize),breaks=seq(0,1000,by=binwidth),plot = F)
      hrc=hist(A2D(subset(RECs,Genus_Code==uG[i_G])$EndingSize),breaks=seq(0,1000,by=binwidth),plot = F)
      PropRec_v=(hrc$counts[1:(N+1)]/hcl$counts[1:(N+1)])
      out_i=which(PropRec_g$Genus_Code==uG[i_G]&PropRec_g$BinWidth==binwidth)
      PropRec_g$PropRec[out_i]=mean(PropRec_v[1:(N+1)],na.rm=T)
      # plot(hcl$breaks[1:(N+1)],PropRec_v,type="b",ylim=c(0,1),main=uG[i_G])
      # abline(h=mean(propREC,na.rm=T))}
    }
  }
}
PropRec_g

detach(package:plyr)
PropRecMn=PropRec_g %>% 
  group_by(Genus_Code) %>% 
  summarise(meanPropRec=mean(PropRec))
#Final 0-5 cm diam. Proportion Recruits in the Juvenile size classes (will use for 'pro-rating')
PropRecMn



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

#Area Surveyed for each SIG
Asurv=read.csv("./data/AreaSurveyed_N_Circrats.csv"); names(Asurv)[1]="Site"
Asurv$Site[Asurv$Site=="OAH_XXX_022"]="OAH_OCC_005"
Asurv$EndingDate=mdy(Asurv$Date)
Asurv_l=Asurv[,c("Site","EndingDate","POCS","POSP","MOSP")]%>%
  pivot_longer(cols=c("POSP","MOSP","POCS"),names_to=c("Genus_Code"),values_to=c("Ncircrats")) %>% 
  mutate(Ncircrats=as.numeric(Ncircrats))
Asurv_l$A.Surv.m2 =Asurv_l$Ncircrats*0.5 #numb of circrats * 0.5m2 (size of circrat) to get area surveyed in m2
Asurv_l$A.Surv.cm2 = Asurv_l$Ncircrats*5000 #area surveyed cm2

#Area of Adult Colonies for each SIG
Atax=ColonyLevel %>%
  group_by(SIG,Site,Interval,Genus_Code,StartingDate,EndingDate,Interval_Years) %>% 
  summarise(AdultCoralArea_cm2=sum(StartingSize)) #Area of adult cm2 = sum all corals in each genus for each year

#Link to Nrecruits from the StartingSize summed area: Adults at beginning of interval generate 
RecSFMTib=ColonyLevel %>%
  filter(TransitionTypeSimple %in% c("RECR")) %>% 
  group_by(SEC_NAME,SIG,Site,Interval,Genus_Code,REGION,StartingDate,EndingDate,Interval_Years) %>% 
  summarize(Nrec=length(which(TransitionTypeSimple=="RECR"))) %>% 
  left_join(Atax) %>% 
  left_join(Asurv_l) %>% 
  mutate(
    # numb recruits/area surveyed cm2
    RecSFM_p_Survcm2=Nrec/A.Surv.cm2,
    #Case 1: Assume Site-Level Stock-Recruitment
    #Recruit N per adult area and annual rate
    RecSFM_p_SiteAdcm2=Nrec/AdultCoralArea_cm2,
    RecSFM_p_SiteAdcm2_Yr=RecSFM_p_SiteAdcm2/Interval_Years) # numb recruits/total adult area cm2 annualized

RecSFMDataFrame=as.data.frame(RecSFMTib)

#save(RecSFMDataFrame, file = paste0("data/",Name,"_RecSFMrates.rdata"))

#plot
hist(RecSFMDataFrame$RecSFM_p_SiteAdcm2_Yr, breaks = 20)


#REGIONAL CALCULATIONS BY GENUS AND REGION
RecSFMDataFrame$RecSFM_p_SiteAdcm2_Yr[RecSFMDataFrame$RecSFM_p_SiteAdcm2_Yr == Inf] <-NA
#Calculate median, mean, and lower/upper quantiles for site-level stock recruitment GENUS values
RecSFMSummary=RecSFMDataFrame %>%
  dplyr::select(SEC_NAME,Site,Interval,Genus_Code,REGION,RecSFM_p_SiteAdcm2_Yr)%>%
  group_by(Genus_Code,REGION)%>%
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
rea_sec=read.csv("data/NWHI_MHI_REA_Data_Sector.csv")

#Get Percent Cover by Genus at each Sector
cov1_sec=read.csv("data/NWHI_MHI_Cover_T1_Data_Sector.csv")
cov3_sec=read.csv("data/NWHI_MHI_Cover_T3_Data_Sector.csv")
uSec=na.omit(unique(ColonyLevel$SEC_NAME))
names(cov1_sec)[2:15]=substr(names(cov1_sec)[2:15],6,999)
meta=names(cov1_sec)[2:6]
classes1=names(cov1_sec)[7:15]
SEclasses1=names(cov1_sec)[22:30]

names(cov3_sec)[2:111]=substr(names(cov3_sec)[2:111],6,999)
classes3=names(cov3_sec)[7:111]
SEclasses3=names(cov3_sec)[118:(104+118)]

#select first 15 cols of cov1_sec and pivot CORAL,MA,TURF,CCA,EMA,SC,I,SED,HAL to
#a new column titled CAT and values to column titled COVER
cov1_lse = cov1_sec[,c(meta,SEclasses1)] %>% 
  pivot_longer(cols = all_of(SEclasses1),names_to="CAT",values_to="SE.COVER")
cov1_lse$CAT=substr(cov1_lse$CAT,10,999) #drop pooled SE
cov1_l = cov1_sec[,c(meta,classes1)] %>% #get cover by category & N
  pivot_longer(cols = all_of(classes1),names_to="CAT",values_to="COVER") %>% 
  left_join(cov1_lse,by=c(meta,"CAT"))
cov1_l=cov1_l[,c("REGION","ISLAND","ANALYSIS_SEC","OBS_YEAR","CAT","COVER","SE.COVER","N")]
cov1_l$SD.COVER=cov1_l$SE.COVER*sqrt(cov1_l$N)
cov1_l$VAR.COVER=cov1_l$SD.COVER^2

cov3_lse = cov3_sec[,c(meta,SEclasses3)] %>% 
  pivot_longer(cols = all_of(SEclasses3),names_to="CAT",values_to="SE.COVER")
cov3_lse$CAT=substr(cov3_lse$CAT,10,999)
cov3_l = cov3_sec[,c(meta,classes3)] %>% 
  pivot_longer(cols = all_of(classes3),names_to="CAT",values_to="COVER") %>% 
  left_join(cov3_lse,by=c(meta,"CAT"))
cov3_l=cov3_l[,c("REGION","ISLAND","ANALYSIS_SEC","OBS_YEAR","CAT","COVER","SE.COVER","N")]
cov3_l$SD.COVER=cov3_l$SE.COVER*sqrt(cov3_l$N)
cov3_l$VAR.COVER=cov3_l$SD.COVER^2


MPPcovers=c("MOBR","MOEN","MONE","MOFO","POCS","POEN","POFO","POMA","PONM","POBR") #cover categories for each genus for Tier3 benthic cat
names(MPPcovers)=c(rep("MOSP",4),"POCS",rep("POSP",5))
MMPlu=names(MPPcovers)
names(MMPlu) = MPPcovers

cov_l=rbind(subset(cov1_l,CAT=="CORAL"),subset(cov3_l,CAT%in%MPPcovers)) 
cov_l$GENUS_CODE=MMPlu[cov_l$CAT] #add genus_code column
cov_l$GENUS_CODE[cov_l$CAT=="CORAL"]="SSSS" #if category is coral, assign SSSS genus code
cov_l$SEC_YEAR=paste0(cov_l$ANALYSIS_SEC,"-",cov_l$OBS_YEAR)
cov_l$REGION=factor(cov_l$REGION,levels=c("MHI","NWHI","PRIAs","MARIAN","SAMOA"))
#SUM all within genus variation at each site
keepcols=c("REGION","ISLAND","ANALYSIS_SEC","OBS_YEAR","GENUS_CODE","SEC_YEAR")
cov_sum=cov_l %>%  #cover by genus code and sector
  group_by(ANALYSIS_SEC,OBS_YEAR,GENUS_CODE) %>% 
  summarize(COVER=sum(COVER),VAR.COVER=sum(VAR.COVER),N.COVER=sum(N))
cov_sum$SD.COVER=sqrt(cov_sum$VAR.COVER)
cov_sum$SE.COVER=cov_sum$SD.COVER/sqrt(cov_sum$N.COVER)

cov=left_join(cov_sum,unique(cov_l[,keepcols]),by=c("ANALYSIS_SEC","OBS_YEAR","GENUS_CODE")) #cov_l by distint benth cats. sum all cover by genus
cov=cov[,c(keepcols,"COVER","SE.COVER","N.COVER")] #reorganize columns
#cover is percent cover (0-100%)


# prorated mean juvenile col density for each sector and each year divided by adult cover sector year genus
#get mean, 5 and 95 quantiles for each sector year genus
#for each of those ^^ 3 numbers, want to prorate depending on genus code. 

table(cov$OBS_YEAR,useNA = "always")
table(rea_sec$ANALYSIS_YEAR,useNA = "always")
names(rea_sec)[which(names(rea_sec)=="Sector")]="ANALYSIS_SEC"
names(rea_sec)[which(names(rea_sec)=="n")]="N.JCD"
rea_sec$OBS_YEAR=rea_sec$ANALYSIS_YEAR

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

#Calculate sector level stock recruitment and propagate error
#sector level recruitment for sectors for this study
Jsec_P_Asec= rea_sec %>%
  filter(ANALYSIS_SEC%in%uSec&GENUS_CODE%in%c("SSSS","MOSP","POSP","POCS")) %>% #filter by target taxa and target sectors (uSec is unique sector list from above)
  select(REGION,ISLAND,ANALYSIS_SEC,OBS_YEAR,GENUS_CODE,Mean_JuvColDen_cm2, SE_JuvColDen_cm2,N.JCD)%>%
  left_join(cov,by=c("REGION","ISLAND","ANALYSIS_SEC","GENUS_CODE","OBS_YEAR")) %>% 
  left_join(PropRecMn,by=c("GENUS_CODE"="Genus_Code")) %>%
  group_by(REGION,ISLAND,ANALYSIS_SEC,OBS_YEAR,GENUS_CODE)%>%
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
  filter(REGION%in%c("NWHI","MHI")&GENUS_CODE%in%c("SSSS","MOSP","POSP","POCS")) %>% #uSec is unique sector list from above
  select(REGION,ISLAND,ANALYSIS_SEC,OBS_YEAR,GENUS_CODE,Mean_JuvColDen_cm2, SE_JuvColDen_cm2,N.JCD)%>%
  left_join(cov,by=c("REGION","ISLAND","ANALYSIS_SEC","GENUS_CODE","OBS_YEAR")) %>% 
  left_join(PropRecMn,by=c("GENUS_CODE"="Genus_Code")) %>%
  group_by(REGION,ISLAND,ANALYSIS_SEC,OBS_YEAR,GENUS_CODE)%>%
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

Jsec_P_ALL %>% #filter(ISLAND%in%c("Hawaii","Maui","Oahu","French Frigate","Lisianski","Kure")) %>% 
  ggplot(aes(MN.RecVal))+geom_histogram(binwidth=0.001)+facet_grid(ISLAND~GENUS_CODE)#+xlim(c(0,.12))
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

Jsec_P_ALL %>% #filter(ISLAND%in%c("Hawaii","Maui","Oahu","French Frigate","Lisianski","Kure")) %>% 
  ggplot(aes(MN.RecVal))+
  geom_histogram(binwidth=0.01)+
  geom_vline(aes(xintercept=MD.RecVal_Sec_All),data=RecVal_Sec_Dists,color="darkcyan",lty=1)+
  geom_vline(aes(xintercept=MD.CI95_LO.RecVal_Sec_All),data=RecVal_Sec_Dists,color="darkcyan",lty=3)+
  geom_vline(aes(xintercept=MD.CI95_HI.RecVal_Sec_All),data=RecVal_Sec_Dists,color="darkcyan",lty=3)+
  facet_grid(REGION~GENUS_CODE)+
  scale_x_sqrt()

#plot site-level stock recruitment values and sector-level stock recruitment values
Site=ggplot(RecSFMDataFrame,aes(RecSFM_p_SiteAdcm2_Yr))+
  geom_histogram()+ theme_bw()+
  facet_grid("Genus_Code")+
  #facet_wrap(Genus_Code~REGION)+
  xlim(c(0,.12)) + xlab("Site-level stock recruitment") # "# recruits  / area adult coral"

Sec= Jsec_P_Asec%>%
  filter(GENUS_CODE%in%c("MOSP","POSP","POCS")) %>%
  ggplot(aes(MN.RecVal))+
  geom_histogram()+theme_bw()+
  facet_grid("GENUS_CODE")+
  #facet_wrap(GENUS_CODE~REGION)+
  xlim(c(0,.12))+ xlab("Proportional sector-level stock recruitment") # "Proportional juv. colony density / cover"

Site/Sec


#Sector Stock Rec for REGIONAL MODEL
Regional_SectorStockRec <- as.data.frame(Jsec_P_Asec) 
Regional_SectorStockRec$MN.RecVal[Regional_SectorStockRec$MN.RecVal == Inf] <- NA
Regional_SectorStockRec$SE.RecVal[Regional_SectorStockRec$SE.RecVal == Inf] <- NA

Regional_SectorStockRec <- Regional_SectorStockRec %>%
  select(ANALYSIS_SEC,OBS_YEAR,GENUS_CODE,REGION,MN.RecVal,SE.RecVal)%>%
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
