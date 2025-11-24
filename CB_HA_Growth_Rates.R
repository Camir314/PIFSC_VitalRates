library(tidyverse)
col=read.csv("./Data/ColonyTransitions/FULLDOMAIN_19.22.23_ColonyTransitions_NoDup.csv")
MDASMA=read.csv("./Data/MetaData/VitalRates_LatLong.csv")
MDASMAkeepcols=c("Region","Island","ESD.Site.Name","Latitude","Longitude")
MDHI=read.csv("./Data/MetaData/HA_Site_Metadata.csv")
MDHI=MDHI[,c("Region","Island_Code","SiteName","Latitude","Longitude")]
names(MDHI)=MDASMAkeepcols
MD=rbind(MDASMA[,MDASMAkeepcols],MDHI)
col=left_join(x=col,y=MD,by=c("REGION"="Region","Island"="Island","Site"="ESD.Site.Name"))

col$GROWTH_LEcm=(col$Shape_Area_ENDcm2-col$Shape_Area_STAcm2)/col$Shape_Leng_STAcm
col$GROWTH_LEcm.yr=col$GROWTH_LEcm/col$Interval_Years
table(col$TransitionType)

GS_HI=col %>% filter(REGION%in%c("MHI","NWHI"),TransitionType==c("GROW","SHRINK")) 

HIlat=GS_HI %>% 
  group_by(Island) %>% summarize(Lat_mn=mean(Latitude,na.rm=T)) %>% arrange(Lat_mn)
GS_HI$Island=factor(GS_HI$Island,levels=HIlat$Island)
GS_HI_IGm=GS_HI %>% group_by(Island,GENUS_CODE) %>% summarize(LE.yr_med=median(GROWTH_LEcm.yr,na.rm=T))

GS_HI %>% 
  ggplot(aes(x=GROWTH_LEcm.yr,fill=Island))+
  geom_histogram(aes(y = ..density..),binwidth = .25)+
  geom_density(binwidth = .1)+
  geom_vline(xintercept = 0,lty=4)+
  geom_vline(aes(xintercept = LE.yr_med),color="red",data=GS_HI_IGm)+
  facet_grid(Island~GENUS_CODE,scales="free_y")+
  theme_bw()

GS_HI %>%
  group_by(Island,Site,GENUS_CODE) %>%
  summarize(LAT_mn=mean(Latitude),q85=quantile(GROWTH_LEcm,.85),ci95=1.96*sd(GROWTH_LEcm)/sqrt(length(GROWTH_LEcm))) %>% 
  ggplot(aes(x=LAT_mn,y=q85,ymin=q85-ci95,ymax=q85+ci95,color=Island))+
  geom_point()+
  geom_errorbar()+
  geom_hline(yintercept = 0)+
  facet_grid("GENUS_CODE")+
  theme_bw()

GS_HI %>% 
  group_by(REGION,Island,Site,GENUS_CODE) %>% summarize(mn_LEcm.yr=mean(GROWTH_LEcm.yr,na.rm=T),
                               se_LEcm.yr=sd(GROWTH_LEcm.yr,na.rm=T)/sqrt(length(GROWTH_LEcm.yr)),
                               LAT_mn=mean(Latitude,na.rm=T)) %>% 
  ggplot(aes(x=LAT_mn,y=mn_LEcm.yr,ymin=mn_LEcm.yr-se_LEcm.yr,ymax=mn_LEcm.yr+se_LEcm.yr,color=Island))+
  geom_point()+
  geom_errorbar()+
  geom_hline(yintercept = 0)+
  facet_grid("GENUS_CODE")+
  theme_bw()

library(lme4)
summary()
mm=lmer(GROWTH_LEcm.yr~GENUS_CODE+1|REGION,data=GS_HI)
library(fitdistrplus)

summary(mm)
GS_HI %>% filter(Interval=="16-19") %>% ggplot(aes(x=Latitude,y=GROWTH_LEcm.yr))+
  geom_point(aes(color=Island))+stat_smooth(method="lm")+facet_wrap("GENUS_CODE")
