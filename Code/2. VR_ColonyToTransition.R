
# Library Calls -----------------------------------------------------------
library(tidyverse)
library(lubridate)
library(factoextra)
library(NbClust)

REGIONYR="ASRAMP23"#"MARAMP22"
# Data Load ---------------------------------------------------------------
Col=read.csv(paste0("./CSV files/ColonyLevel/",REGIONYR,"_VitalRates_colonylevel_CLEAN.csv"))
Col$TL_Date=ymd(Col$TL_Date)

# Determine Unique Time Points ---------------------------
#Use YEAR OR DATE_CLUSTER
TP_METHOD="DATE_CLUSTER"
if(TP_METHOD=="YEAR"){
  # It would be easiest to assume a single timepoint per sampling year,
  uY=sort(unique(Col$Year))
  # Make simple lookup
  Year_TP_LU=1:length(uY);names(Year_TP_LU)=uY
  #Assign TimePoint IDs
  Col$TP_ID=Year_TP_LU[as.character(Col$Year)]
  
}else if(TP_METHOD=="DATE_CLUSTER"){
  # but worth checking clustering at Date Level
  uD=sort(unique(Col$TL_Date))
  uDD=dist(uD)
  #Determine best number of TimePoints by clustering
  # K=NbClust(data = uD,diss = uDD,distance = NULL,
  #          min.nc = 2, max.nc = 10, method = "complete")
  # 
  ################################################################################
  ################################################################################
  #Eyeball because Nbclust is failing
  ################################################################################
  ################################################################################
  Kval=3
  #Kval=max(K$Best.partition)
  
  #Check out dendrogram
  hcuDD=hclust(uDD)
  par(mfrow=c(1,1))
  plot(hcuDD,labels=uD);rect.hclust(hcuDD,k=Kval)
  
  #write Date to Timepoint number Look up table
  Date_TP_df=data.frame(ID0=cutree(hcuDD,k=Kval),Date=uD)
  #ensure TP # are chronological
  ChronID=Date_TP_df %>% group_by(ID0) %>% summarize(mnDate=mean(Date)) %>% arrange(mnDate) 
  ChronID$IDc=1:nrow(ChronID)
  
  #Build Chronological Lookup
  Date_TP_df=left_join(Date_TP_df,ChronID[,c("ID0","IDc")],by="ID0")
  Date_TP_LU=Date_TP_df$IDc;names(Date_TP_LU)=Date_TP_df$Date
  
  #Assign TimePoint IDs
  Col$TP_ID=Date_TP_LU[as.character((Col$TL_Date))]
}else{
  stop("No valid Time Point selection method.")
}
#table(Col$TP_ID,Col$Site)


# Two-Stage Transition Loop - Rbind Aggregation ---------------------------
uTP=sort(unique(Col$TP_ID))
#Loop through one less than the total number of TPs
#site data
SiteData=unique(Col[,c("Island","Site","TP_ID","Year","TL_Date")]) %>% arrange(Site,TP_ID)

#joincols
alljoincols=c("Island","Site","Genus","Site_Genet","TP_ID","Year","TL_Date","TL_Area","TL_Perim","nPatches")
byjoincols=c("Island","Site","Genus","Site_Genet")

ColonyTransitions=NULL
for(i_tp in 1:(length(uTP)-1)){
  TP_0=i_tp;TP_1=i_tp+1
  
  #For the data across these two time points, only include sites sampled in both years
  TP_0sites=SiteData %>% filter(TP_ID == TP_0) %>% select(Site) %>% distinct()
  TP_1sites=SiteData %>% filter(TP_ID == TP_1) %>% select(Site) %>% distinct()
  SampledSites=intersect(TP_0sites,TP_1sites)
  #Subet colonies from sites that span the interval
  Col_tp2=Col %>% 
    filter(TP_ID %in% c(TP_0,TP_1)) %>% 
    filter(Site %in% SampledSites$Site)
  
  #Count Colonies
  tSG=table(Col_tp2$Site_Genet)  #check: table(tSG)
  
  #Colonies with GROW Transitions
  GTr=names(tSG[which(tSG==2)])
  GTrCol_0=Col_tp2 %>% 
    filter(Site_Genet %in% GTr) %>% 
    filter(TP_ID == TP_0) %>% 
    select(all_of(alljoincols))
  GTrCol_1=Col_tp2 %>% 
    filter(Site_Genet %in% GTr) %>% 
    filter(TP_ID == TP_1) %>% 
    select(all_of(alljoincols))
  GTrCol=left_join(GTrCol_0,GTrCol_1,by=byjoincols,suffix=c("_STA","_END"))
  GTrCol$TransitionTypeSimple="GROW"
  
  #Colonies with RECR/MORT Transitions
  RMTr=names(tSG[which(tSG==1)])
  
  #extract data for RECR Trans
  RTrCol_1=Col_tp2 %>% 
    filter(Site_Genet %in% RMTr) %>% 
    filter(TP_ID == TP_1)  %>% 
    select(all_of(alljoincols))
  #Build Fake TP_0 Dataframe
  RTrCol_0=RTrCol_1[,byjoincols]
  RTrCol_0$TP_ID=TP_0
  RTrCol_0=left_join(RTrCol_0,SiteData,by=c("Island","Site","TP_ID"))
  RTrCol_0$TL_Area=0
  RTrCol_0$TL_Perim=0
  RTrCol_0$nPatches=0
  RTrCol=left_join(RTrCol_0,RTrCol_1,by=byjoincols,suffix=c("_STA","_END"))
  RTrCol$TransitionTypeSimple="RECR"
  
  #extract data for MORT Trans
  MTrCol_0=Col_tp2 %>% 
    filter(Site_Genet %in% RMTr) %>% 
    filter(TP_ID == TP_0)  %>% 
    select(all_of(alljoincols))
  #Build Fake TP_1 Dataframe
  MTrCol_1=MTrCol_0[,byjoincols]
  MTrCol_1$TP_ID=TP_1
  MTrCol_1=left_join(MTrCol_1,SiteData,by=c("Island","Site","TP_ID"))
  MTrCol_1$TL_Area=0
  MTrCol_1$TL_Perim=0
  MTrCol_1$nPatches=0
  MTrCol=left_join(MTrCol_0,MTrCol_1,by=byjoincols,suffix=c("_STA","_END"))
  MTrCol$TransitionTypeSimple="MORT"
  
  TP2_Trans=rbind(GTrCol,RTrCol,MTrCol)
  ColonyTransitions=rbind(ColonyTransitions,TP2_Trans)
}

#Build Out Colony Transition Data
ColonyTransitions=ColonyTransitions %>%
  mutate(
    #Interval Name
    Interval=paste0(substr(Year_STA,3,4),"_",substr(Year_END,3,4)),
    #Transition and Annualization
    l10_Area_STA=log10(TL_Area_STA),
    l10_Area_END=log10(TL_Area_END),
    l10TransitionMagnitude=l10_Area_END-l10_Area_STA,
    Interval_Years=as.numeric(difftime(TL_Date_END,TL_Date_STA,units = "days"))/365.25,
    l10_Area_ENDann=l10_Area_STA+(l10TransitionMagnitude/Interval_Years),
    #Mortality and Survival
    Mortality=if_else(TransitionTypeSimple=="MORT",1,0),
    Survival=1-Mortality
)

# Transition Data Out -----------------------------------------------------
write.csv(x = ColonyTransitions,file = paste0("./CSV files/ColonyTransitions/",REGIONYR,"_ColonyTransitions.csv"))

