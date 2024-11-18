library(tidyverse)
COV=read.csv("T:/Benthic/Data/REA Coral Demography & Cover/Summary Data/Island/BenthicCover_2010-2019_Tier1_ISLAND.csv")
REA=read.csv("T:/Benthic/Data/REA Coral Demography & Cover/Summary Data/Island/BenthicREA_islanddata_TAXONCODE_updated.csv")
LU=read.csv("T:/Benthic/Data/Lookup Tables/2013-23_Taxa_MASTER.csv")

FFSc=COV %>% filter(Mean.ISLAND=="French Frigate")
FFSr_all=REA %>% filter(Island=="French Frigate")
FFSr_SSSS=FFSr_all %>% filter(TAXONCODE=="SSSS")

ggplot(FFSr_SSSS,
       aes(x=ANALYSIS_YEAR,y=Mean_AdColDen,ymin=Mean_AdColDen-SE_AdColDen,ymax=Mean_AdColDen+SE_AdColDen ))+
  geom_col(color="black",fill="lightblue")+
  geom_errorbar(width=.5)+
  theme_bw()+ggtitle("Lalo: Adult Colony Density (All Corals)")+
  xlab("Analysis Year")+
  ylab("Colony Density (col/m2)")


ggplot(FFSc,aes(x=Mean.ANALYSIS_YEAR,y=Mean.CORAL,ymin=Mean.CORAL-PooledSE.CORAL,ymax=Mean.CORAL+PooledSE.CORAL ))+
  geom_col(color="black",fill="lightblue")+
  geom_errorbar(width=.5)+
  theme_bw()+ggtitle("Lalo: Hard Coral Percent Cover (All Corals)")+
  xlab("Analysis Year")+
  ylab("Hard Coral (% Cover)")


TaxAll=FFSr_all %>% filter(Mean_AdColDen>0.025,ANALYSIS_YEAR=="2016") %>% 
  group_by(TAXONCODE) %>%
  summarize(ACD=(Mean_AdColDen)) %>% 
  #  filter(TAXONCODE!="SSSS") %>% 
  arrange(desc(ACD)) %>% 
  mutate(TaxonName=LU$TAXON_NAME[match(TAXONCODE,LU$SPCODE)]) 
TaxAll=TaxAll[,c(3,1,2)]
TaxAll$TaxonName[1]="All Corals"
TaxAll$TaxonName=factor(TaxAll$TaxonName,levels=TaxAll$TaxonName)
#FFStaxa=LU$TAXON_NAME[match(TaxAll$TAXONCODE,LU$SPCODE)]
#write.csv(TaxAll,"T:/Benthic/Data/Lookup Tables/FFStaxa.csv")

ggplot(TaxAll,aes(x=TaxonName,y=ACD))+
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  ggtitle("Lalo: Adult Colony Density (By Taxon)")+
  xlab("Taxon Name")+
  ylab("Colony Density (col/m2)")

ggplot(REA_S,aes(x=AdColDen))+
  geom_histogram()+
  scale_x_log10()+
  facet_wrap("ISLAND",scales="free_y")
min(REA_S$AdColDen[which(REA_S$AdColDen>0)],na.rm=T)
