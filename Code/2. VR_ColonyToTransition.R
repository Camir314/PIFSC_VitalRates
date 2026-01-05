
# Library Calls -----------------------------------------------------------
library(tidyverse)
library(lubridate)
library(factoextra)
library(NbClust)

REGIONYR=c("MARAMP22","ASRAMP23","HARAMP24")# ##
for (iRY in REGIONYR){
  # Data Load ---------------------------------------------------------------
  Col=read.csv(paste0("./Data/ColonyLevel/",iRY,"_VitalRates_colonylevel_CLEAN.csv"))
  Col$TL_Date=ymd(Col$TL_Date)
  Col=Col %>% arrange(Site_Genet,TimePt) %>% filter(TL_Class != "Dummy")
  if (iRY%in%c("MARAMP22","ASRAMP23")){Col=Col %>% rename(Surface_Area=SArea)}
    
  # Determine Unique Time Points ---------------------------
  #Use YEAR OR DATE_CLUSTER
  TP_METHOD="YEAR"#"DATE_CLUSTER"
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
    
    #write Date to Time point number Look up table
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
  uTP=sort(unique(Col$TimePt))
  #Loop through one less than the total number of TPs
  #site data
  SiteData=unique(Col[,c("Island","Site","TimePt","TP_ID","Year","TL_Date")]) %>% arrange(Site,TP_ID)
  
  #joincols
  alljoincols=c("Island","Site","Genus","Site_Genet","TimePt","TP_ID","Year","TL_Date","Surface_Area","Shape_Area","Shape_Leng","nPatches") #switch to shape_area and shape_perim 
  byjoincols=c("Island","Site","Genus","Site_Genet")
  
  #Search for raw dups and drop (should be none)
  Dup0_i=Col %>% duplicated(nmax=nrow(Col)) %>% which()
  Col=Col %>% filter(!(Genet_full %in% Col$Genet_full[Dup0_i]))
  
  #Search duplicates that differ only by quadrat (differing by Quadrat #)
  DupQ_i=Col %>% select(-Quadrat) %>% duplicated(nmax=nrow(Col)) %>% which()
  Col %>% filter(Genet_full%in% Col$Genet_full[DupQ_i])
  #Drops any Quadrat Dups (losing -99)
  Col=Col %>% filter(!(Genet_full %in% Col$Genet_full[DupQ_i]&Quadrat==-99))
  
  #Search duplicates that differ only by Species, if genus is same, keep one, if genus different (2) drop em both
  DupS_i=Col %>% select(-TL_Class) %>% duplicated(nmax=nrow(Col)) %>% which()
  #Dup_i=Col %>% select(Shape_Leng,Shape_Area,Surface_Area) %>% duplicated(nmax=nrow(Col)) %>% which()
  Col %>% filter(Genet_full%in% Col$Genet_full[DupS_i])
  
  #First change TL_Class to Genus
  Col=Col %>% mutate(TL_Class=ifelse((Genet_full %in% Col$Genet_full[DupS_i]),
                                     Genus,TL_Class))
  #Now that TL_Class is same as Genus for duplicated rows, drop any overall dups
  Dup02_i=Col %>% duplicated(nmax=nrow(Col)) %>% which()
  if(length(Dup02_i)>0){Col=Col[-Dup02_i,]}
  
  #Now Get genus Dups, drop both entries
  DupG_i=Col %>% select(-Genus,-TL_Class) %>% duplicated(nmax=nrow(Col)) %>% which()
  Col %>% filter(Genet_full%in% Col$Genet_full[DupG_i])
  Col=Col %>% filter(!(Genet_full%in% Col$Genet_full[DupG_i]))
  
  #Check Taxon ID for changes
  SplitTax_sg=Col %>% group_by(Site_Genet) %>% summarize(Ntax=length(unique(Genus))) %>% filter(Ntax>1) %>% pull(Site_Genet)
  Col=Col %>% filter(!Site_Genet%in%SplitTax_sg)
  
  #This code is structured to loop timepoints and quickly include all single-timestep
  #Transitions. In rare cases we may mistake a missed timepoint for a MORT then RECR
  #I'll trouble shoot that at the end.
  ColonyTransitions=NULL
  for(i_tp in 1:(length(uTP)-1)){
    TP_0=uTP[i_tp];    TP_1=TP_0+1;
    # TP_0=2; TP_1=3
    
    #For the data across these two time points, only include sites sampled in both years
    TP_0sites=SiteData %>% filter(TimePt == TP_0) %>% select(Site) %>% distinct()
    TP_1sites=SiteData %>% filter(TimePt == TP_1) %>% select(Site) %>% distinct()
    SampledSites=intersect(TP_0sites,TP_1sites)
    
    #If no site share those time point, advance tp_1 by one...
    if(nrow(SampledSites)==0){next}
    
    #Subet colonies from sites that span the interval
    Col_tp2=Col %>% 
      filter(TimePt %in% c(TP_0,TP_1)) %>% 
      filter(Site %in% SampledSites$Site)
    
    #Count Colonies
    tSG=table(Col_tp2$Site_Genet)  #check: table(tSG)
    #table(tSG)
    #Colonies with GROW Transitions
    GTr=names(tSG[which(tSG==2)])
    # #Drop Colonies that have intermediate timepoint values
    # if(TP_1!=(TP_0+1)){
    #   DropMultiSpanCol_GR = Col %>% filter(Site_Genet %in% GTr,TP_ID %in% (TP_0+1):(TP_1-1)) %>% pull(Site_Genet)
    # }else{
    #   DropMultiSpanCol=NULL
    # }
    
    #Get TP_0 data
    GTrCol_0=Col_tp2 %>% 
      filter(Site_Genet %in% GTr) %>% #,!Site_Genet %in% DropMultiSpanCol) %>% 
      filter(TimePt == TP_0) %>% 
      select(all_of(alljoincols))
    GTrCol_1=Col_tp2 %>% 
      filter(Site_Genet %in% GTr) %>% #,!Site_Genet %in% DropMultiSpanCol) %>% 
      filter(TimePt == TP_1) %>% 
      select(all_of(alljoincols))
    GTrCol=left_join(GTrCol_0,GTrCol_1,by=byjoincols,suffix=c("_STA","_END"))
    if(nrow(GTrCol)==0){GTrCol$TransitionTypeSimple=NULL}else{GTrCol$TransitionTypeSimple="GROW"}
    
    #Colonies with RECR/MORT Transitions
    RMTr=names(tSG[which(tSG==1)])
    
    #extract data for RECR Trans
    RTrCol_1=Col_tp2 %>% 
      filter(Site_Genet %in% RMTr) %>% #,!Site_Genet %in% DropMultiSpanCol) %>% 
      filter(TimePt == TP_1)  %>% 
      select(all_of(alljoincols))
    #Build Fake TP_0 Dataframe
    RTrCol_0=RTrCol_1[,byjoincols]
    RTrCol_0$TimePt=TP_0
    RTrCol_0=left_join(RTrCol_0,SiteData,by=c("Island","Site","TimePt"))
    RTrCol_0$Shape_Area=0
    RTrCol_0$Shape_Leng=0
    RTrCol_0$Surface_Area=0
    RTrCol_0$nPatches=0
    RTrCol=left_join(RTrCol_0,RTrCol_1,by=byjoincols,suffix=c("_STA","_END"))
    if(nrow(RTrCol)==0){RTrCol$TransitionTypeSimple=NULL}else{RTrCol$TransitionTypeSimple="RECR"}
    
    #extract data for MORT Trans
    MTrCol_0=Col_tp2 %>% 
      filter(Site_Genet %in% RMTr) %>% #,!Site_Genet %in% DropMultiSpanCol) %>% 
      filter(TimePt == TP_0)  %>% 
      select(all_of(alljoincols))
    #Build Fake TP_1 Dataframe
    MTrCol_1=MTrCol_0[,byjoincols]
    MTrCol_1$TimePt=TP_1
    MTrCol_1=left_join(MTrCol_1,SiteData,by=c("Island","Site","TimePt"))
    MTrCol_1$Shape_Area=0
    MTrCol_1$Shape_Leng=0
    MTrCol_1$Surface_Area=0
    MTrCol_1$nPatches=0
    MTrCol=left_join(MTrCol_0,MTrCol_1,by=byjoincols,suffix=c("_STA","_END"))
    if(nrow(MTrCol)==0){MTrCol$TransitionTypeSimple=NULL}else{MTrCol$TransitionTypeSimple="MORT"}
    
    TP2_Trans=rbind(GTrCol,RTrCol,MTrCol) 
    ColonyTransitions=rbind(ColonyTransitions,TP2_Trans)
  }
  
  #Find any Transitions where MORT isn't the last for a Site_Genet, or RECR isn't the first
  BrokenTrans=ColonyTransitions %>%
    group_by(Site_Genet) %>%
    mutate(firstTP=min(TimePt_STA),lastTP=max(TimePt_STA),
           RTP=ifelse(TransitionTypeSimple==c("RECR"),TimePt_STA!=firstTP,FALSE),
           MTP=ifelse(TransitionTypeSimple==c("MORT"),TimePt_STA!=lastTP,FALSE)) %>%
    filter(RTP==TRUE|MTP==TRUE) #%>% View()
  
  #Remove 'em from CT
  CTminus=ColonyTransitions %>%
    group_by(Site_Genet) %>%
    mutate(firstTP=min(TimePt_STA),lastTP=max(TimePt_STA),
           RTP=ifelse(TransitionTypeSimple==c("RECR"),TimePt_STA!=firstTP,FALSE),
           MTP=ifelse(TransitionTypeSimple==c("MORT"),TimePt_STA!=lastTP,FALSE)) %>%
    filter(RTP==FALSE,MTP==FALSE)# %>% View()
  dim(ColonyTransitions);dim(BrokenTrans);dim(CTminus);
  
  #Fix em
  IDcol=c("Island", "Site"  ,   "Genus", "Site_Genet")
  STACol=grep(pattern = "_STA",names(ColonyTransitions),value = TRUE)
  ENDCol=grep(pattern = "_END",names(ColonyTransitions),value = TRUE)
  BT_mort=BrokenTrans %>% filter(TransitionTypeSimple=="MORT")
  BT_recr=BrokenTrans %>% filter(TransitionTypeSimple=="RECR")
  UnBrokenTrans=BrokenTrans %>% select(IDcol) %>% distinct() %>% 
    left_join(BT_mort[c(IDcol,STACol)],by=all_of(IDcol)) %>% 
    left_join(BT_recr[c(IDcol,ENDCol)],by=all_of(IDcol)) 
  UnBrokenTrans$TransitionTypeSimple="GROW"
  
  ColonyTransitions.=rbind(CTminus,UnBrokenTrans) %>%
    select(-firstTP,-lastTP,-RTP,-MTP) %>% 
    arrange(Site_Genet,TimePt_STA)
  dim(ColonyTransitions);dim(ColonyTransitions.)
  ColonyTransitions=ColonyTransitions.;rm(list="ColonyTransitions.")
  
  #Build Out Colony Transition Data
  ColonyTransitions=ColonyTransitions %>%
    mutate(
      #Interval Name
      Interval=paste0(substr(Year_STA,3,4),"_",substr(Year_END,3,4)),
      #Transition and Annualization
      l10_Area_STA=log10(Shape_Area_STA),
      l10_Area_END=log10(Shape_Area_END),
      l10TransitionMagnitude=l10_Area_END-l10_Area_STA,
      l10_Surface_Area_STA=log10(Surface_Area_STA),
      l10_Surface_Area_END=log10(Surface_Area_END),
      l10STransitionMagnitude=l10_Surface_Area_END-l10_Surface_Area_STA,
      Interval_Years=as.numeric(difftime(TL_Date_END,TL_Date_STA,units = "days"))/365.25,
      l10_Area_ENDann=l10_Area_STA+(l10TransitionMagnitude/Interval_Years),
      l10_Surface_Area_ENDann=l10_Surface_Area_STA+(l10STransitionMagnitude/Interval_Years),
      #Mortality and Survival
      Mortality=if_else(TransitionTypeSimple=="MORT",1,0),
      Survival=1-Mortality
    )
  #Check Intervals
  sort(unique(ColonyTransitions$Interval))
  plot((table(ColonyTransitions$Interval_Years)))
  ColonyTransitions %>% filter(Interval_Years<1) %>% as.data.frame()
  
  # Transition Data Out -----------------------------------------------------
  print(paste0("Writing Region Year: ",iRY))
  write.csv(x = ColonyTransitions,file = paste0("./Data/ColonyTransitions/",iRY,"_ColonyTransitions.csv"))
}


# Compile Transition Data Out -----------------------------------------------------
ColTransPath="./Data/ColonyTransitions/"
RegionalFiles=sort(list.files(path = ColTransPath,pattern = "_ColonyTransitions.csv"))
MA=read.csv(paste0(ColTransPath,RegionalFiles[3]))
AS=read.csv(paste0(ColTransPath,RegionalFiles[1]))
HA=read.csv(paste0(ColTransPath,RegionalFiles[2]))
MA=MA %>% select(-X,-TimePt_STA,-TimePt_END)  
AS=AS %>% select(-X,-TimePt_STA,-TimePt_END) 
HA=HA %>% select(-X,-TP_ID_STA,-TP_ID_END) %>%
  rename(TP_ID_STA=TimePt_STA,TP_ID_END=TimePt_END)

names(MA)
names(AS)
names(HA)
setdiff(names(MA),names(AS))
setdiff(names(AS),names(MA))
setdiff(names(MA),names(HA))
setdiff(names(HA),names(MA))

MAASHA=rbind(MA,AS,HA)
uI=sort(unique(MAASHA$Island))
RegLU=c("MARIAN","PRIA","NWHI","MARIAN","MHI","PRIA","MHI","MARIAN","MHI","SAMOA","MARIAN","NWHI","SAMOA","MARIAN","SAMOA","SAMOA")
names(RegLU)=uI
MAASHA$REGION=RegLU[MAASHA$Island]


#
write.csv(x = MAASHA,file = paste0("./Data/ColonyTransitions/MAASHA_22-24_AllColonyTransitions.csv"),row.names = FALSE)


dim(MAASHA)
length(unique(MAASHA$Site_Genet))
table(MAASHA$Site,MAASHA$Interval,MAASHA$Genus,MAASHA$TransitionTypeSimple)

uSG=unique(ColonyTransitions$Site_Genet)
thisgenet=sample(uSG,1);thisgenet
thisgenet="OCC-MAI-017_490"
ColonyTransitions %>% filter(Site_Genet==thisgenet) %>% as.data.frame()


ggg=ColonyTransitions %>% 
  group_by(Site_Genet) %>% 
  summarize(G.="GROW"%in%TransitionTypeSimple,M.="MORT"%in%TransitionTypeSimple,R.="RECR"%in%TransitionTypeSimple,
            NTrans=length(TransitionTypeSimple),RGM.=R.&G.&M.) %>% arrange(desc(RGM.),desc(NTrans))

hist(ggg$NTrans)
