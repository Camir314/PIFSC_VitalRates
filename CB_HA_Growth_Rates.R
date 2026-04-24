rm(list=ls())
library(patchwork)
library(tidyverse)
LE_vol=function(SA1,PA1,P1,SA2,PA2,P2){
  return(((PA2*((SA2-PA2)/P2))-(PA1*((SA1-PA1)/P1)))/SA1)
}

# col=read.csv("./Data/ColonyTransitions/FULLDOMAIN_19.22.23_ColonyTransitions_NoDup.csv")
MASHA=read.csv("./Data/ColonyTransitions/MAASHA19-24_ColonyTransitions.csv")
col=MASHA %>% filter(Island %in% c("FFS","HAW","KUR","LIS","MAI","OAH","PHR"))
col$TransitionTypeSimple[col$TransitionTypeSimple=="GROWTH"]="GROW"
col$TransitionTypeSimple[col$TransitionTypeSimple=="SHRINK"]="GROW"

MDASMA=read.csv("./Data/MetaData/VitalRates_LatLong.csv")
MDASMAkeepcols=c("Region","Island","Site","Latitude","Longitude")

MDHI=read.csv("./Data/MetaData/HA_Site_Metadata.csv")
MDHI=MDHI[,c("Region","Island_Code","SiteName","Latitude","Longitude")]
names(MDHI)=MDASMAkeepcols
MD=rbind(MDASMA[,MDASMAkeepcols],MDHI[,MDASMAkeepcols]) %>% arrange(Site) %>% distinct()
col=left_join(x=col,y=MD,by=c("Island"="Island","Site"))

length(which(col$Surface_Area_STA==0))
dim(col)
gHI=col %>% 
  filter(TransitionTypeSimple=="GROW",Region%in%c("MHI","NWHI"),
         Surface_Area_STA!=0,Surface_Area_END!=0,
         !is.na(Surface_Area_STA),!is.infinite(Surface_Area_STA),
         !is.na(Shape_Area_STA),!is.infinite(Shape_Area_STA),
         !is.na(Surface_Area_END),!is.infinite(Surface_Area_END),
         !is.na(Shape_Area_END),!is.infinite(Shape_Area_END)) %>%  
  mutate(Perimeter_STAcm=100*Shape_Leng_STA,
         Planar_Area_STAcm2=10000*Shape_Area_STA,
         Surface_Area_STAcm2=10000*Surface_Area_STA,
         Ave_H_STAcm=(Surface_Area_STAcm2-Planar_Area_STAcm2)/Perimeter_STAcm,
         Vol_STAcm3=Planar_Area_STAcm2*Ave_H_STAcm,
         Perimeter_ENDcm=100*Shape_Leng_END,
         Planar_Area_ENDcm2=10000*Shape_Area_END,
         Surface_Area_ENDcm2=10000*Surface_Area_END,
         Ave_H_ENDcm=(Surface_Area_ENDcm2-Planar_Area_ENDcm2)/Perimeter_ENDcm,
         Vol_ENDcm3=Planar_Area_ENDcm2*Ave_H_ENDcm,
         d_Pcm.yr=(Perimeter_ENDcm-Perimeter_STAcm)/Interval_Years,
         d_PAcm2.yr=(Planar_Area_ENDcm2-Planar_Area_STAcm2)/Interval_Years,
         d_SAcm2.yr=(Surface_Area_ENDcm2-Surface_Area_STAcm2)/Interval_Years,
         dL_Pcm.yr=(log(Perimeter_ENDcm)-log(Perimeter_STAcm))/Interval_Years,
         dL_PAcm2.yr=(log(Planar_Area_ENDcm2)-log(Planar_Area_STAcm2))/Interval_Years,
         dL_SAcm2.yr=(log(Surface_Area_ENDcm2)-log(Surface_Area_STAcm2))/Interval_Years,
         d_Hcm.yr=(Ave_H_ENDcm-Ave_H_STAcm)/Interval_Years,
         d_Vcm3.yr=(Vol_ENDcm3-Vol_STAcm3)/Interval_Years,
         d_nPatches=nPatches_END-nPatches_STA,
         LinExt_Planar_cm=(Planar_Area_ENDcm2-Planar_Area_STAcm2)/Perimeter_STAcm,
         LinExt_Planar_cm.yr=LinExt_Planar_cm/Interval_Years,
         LinExt_Vol_cm=LE_vol(SA1 = Surface_Area_STAcm2,SA2 = Surface_Area_ENDcm2,
                              PA1 = Planar_Area_STAcm2,PA2=Planar_Area_ENDcm2,
                              P1 = Perimeter_STAcm,P2 = Perimeter_ENDcm),
         LinExt_Vol_cm_CHECK=(Vol_ENDcm3-Vol_STAcm3)/Surface_Area_STAcm2,
         LinExt_Vol_cm.yr=LinExt_Vol_cm/Interval_Years,
         LinExt_Vol_cm.yr_CHECK=d_Vcm3.yr/Surface_Area_STAcm2
  ) %>% 
  filter(!is.na(LinExt_Vol_cm.yr),!is.infinite(LinExt_Vol_cm.yr)) %>%
  dplyr::select(Island,Site,Genus_Code,Site_Genet,Latitude,Longitude,
                Year_STA,Date_STA,Year_END,Date_END,Interval,Interval_Years,
                Perimeter_STAcm,Planar_Area_STAcm2,Surface_Area_STAcm2,nPatches_STA,
                Perimeter_ENDcm,Planar_Area_ENDcm2,Surface_Area_ENDcm2,nPatches_END,
                d_Pcm.yr,d_PAcm2.yr,d_SAcm2.yr,dL_Pcm.yr,dL_PAcm2.yr,dL_SAcm2.yr,
                d_Hcm.yr,d_Vcm3.yr,d_nPatches,
                LinExt_Planar_cm,LinExt_Vol_cm,
                LinExt_Planar_cm.yr,LinExt_Vol_cm.yr)

table(gHI$Genus_Code,gHI$Island)

LEvVSp=gHI %>% ggplot(aes(x=LinExt_Planar_cm.yr,y=LinExt_Vol_cm.yr,color=Genus_Code))+
  geom_point()+scale_y_log10()+scale_x_log10()+geom_abline()+theme_bw();LEvVSp

qUPv=quantile(gHI$LinExt_Vol_cm.yr,.975)
qUPp=quantile(gHI$LinExt_Planar_cm.yr,.975)
qDNv=quantile(gHI$LinExt_Vol_cm.yr,.025)
qDNp=quantile(gHI$LinExt_Planar_cm.yr,.025)
PAvLEv=gHI %>%
  filter(LinExt_Vol_cm.yr<qUPv) %>%
  filter(LinExt_Vol_cm.yr>qDNv) %>%
  #  filter(d_nPatches!=0) %>%
  ggplot(aes(x=Planar_Area_STAcm2,y=LinExt_Vol_cm.yr,color=d_nPatches))+
  geom_text(size=3,alpha=1,aes(label=(d_nPatches)))+
  #geom_point(size=3,alpha=.5,aes(size=abs(d_nPatches)))+
  scale_x_log10()+geom_hline(yintercept=0)+
  scale_color_viridis_c(limits=c(-5,5),oob=squish)+
  geom_hline(yintercept=qDNv,lty=2)+
  geom_hline(yintercept=qUPv,lty=2)+theme_bw()+
  facet_wrap("Genus_Code");PAvLEv

PAvLEp=gHI %>%
  filter(LinExt_Planar_cm.yr<qUPp) %>%
  filter(LinExt_Planar_cm.yr>qDNp) %>%
  #  filter(d_nPatches!=0) %>%
  ggplot(aes(x=Planar_Area_STAcm2,y=LinExt_Planar_cm.yr,color=d_nPatches))+
  geom_text(size=3,alpha=1,aes(label=(d_nPatches)))+
  #geom_point(size=3,alpha=.5,aes(size=abs(d_nPatches)))+
  scale_x_log10()+geom_hline(yintercept=0)+
  scale_color_viridis_c(limits=c(-5,5),oob=squish)+
  geom_hline(yintercept=qUPp,lty=2)+
  geom_hline(yintercept=qDNp,lty=2)+
  theme_bw()+
  facet_wrap("Genus_Code");PAvLEp

PAvLEp/PAvLEv


HIlat=gHI %>%
  group_by(Island) %>% summarize(Lat_mn=mean(Latitude,na.rm=T)) %>% arrange(Lat_mn)
gHI$Island=factor(gHI$Island,levels=HIlat$Island)

gHI.=gHI %>% filter(LinExt_Vol_cm.yr<qUPv,LinExt_Vol_cm.yr>qDNv)

gHI_IGm=gHI. %>% group_by(Island,Genus_Code) %>%
  summarize(LE.yr_med=median(LinExt_Vol_cm.yr,na.rm=T))


gHI. %>% 
  ggplot(aes(x=LinExt_Vol_cm.yr,fill=Island))+
  geom_histogram(aes(y = ..density..),binwidth = .25)+
  xlim(c(1.1*qDNv,1.1*qUPv))+
  geom_density(binwidth = .025)+
  geom_vline(xintercept = c(qDNv,qUPv),lty=4)+
  geom_vline(xintercept = c(0),lty=1)+
  geom_vline(aes(xintercept = LE.yr_med),color="red",data=gHI_IGm)+
  facet_grid(Island~Genus_Code,scales="free_y")+
  theme_bw()

gHI. %>%
  group_by(Island,Site,Genus_Code) %>%
  summarize(LAT_mn=mean(Latitude),
            N=length(LinExt_Vol_cm.yr),
            q85=quantile(LinExt_Vol_cm.yr,.5),
            ci95=qt(.975,nrow(gHI.))*sd(LinExt_Vol_cm.yr)/sqrt(length(LinExt_Vol_cm.yr))) %>% 
  ggplot(aes(x=LAT_mn,y=q85,ymin=q85-ci95,ymax=q85+ci95,color=Island))+
  geom_point(aes(size=N))+
  geom_errorbar()+
  geom_hline(yintercept = 0)+
  facet_grid("Genus_Code",scales="free_y")+
  theme_bw()

gHI. %>% 
  group_by(Island,Site,Genus_Code) %>% 
  summarize(mn_LEcm.yr=mean(LinExt_Vol_cm.yr,na.rm=T),
            se_LEcm.yr=sd(LinExt_Vol_cm.yr,na.rm=T)/sqrt(length(LinExt_Vol_cm.yr)),
            LAT_mn=mean(Latitude,na.rm=T)) %>% 
  ggplot(aes(x=LAT_mn,y=mn_LEcm.yr,ymin=mn_LEcm.yr-se_LEcm.yr,ymax=mn_LEcm.yr+se_LEcm.yr,color=Island))+
  geom_point()+
  geom_errorbar()+
  geom_hline(yintercept = 0)+
  facet_grid("Genus_Code",scales="free_y")+
  theme_bw()

library(lme4)
summary()
mm=lmer(LinExt_Vol_cm.yr~Genus_Code+1|Island,data=gHI.)

library(fitdistrplus)

summary(mm)
gHI. %>% filter(Interval=="16_19") %>% ggplot(aes(x=Latitude,y=LinExt_Vol_cm.yr))+
  geom_point(aes(color=Island))+stat_smooth(method="lm")+facet_wrap("Genus_Code")

dim(gHI.)
table(gHI.$Island,gHI.$Genus_Code)

write.csv(gHI.,"./Data/Hawaii_LE_Rates_For_Carb_Budg_Centralq95.csv",row.names = F)

#########################
#Experiment with Growth LE metrics

LE_vol=function(SA1,PA1,P1,SA2,PA2,P2){
  return(((PA2*((SA2-PA2)/P2))-(PA1*((SA1-PA1)/P1)))/SA1)
}

vs=read.csv("./Data/VolumeScenarios/VolScen.csv")
vs
N=choose(n = nrow(vs),k=2)
tr=data.frame(scen=rep(NA,N))#,delPA=NA,delSA=NA,delP=NA,delV=NA,delH=NA,PlanarLE=NA,up=NA,delRug=NA,sap=NA,VolumeLE_proxy=NA,VolumeLE=NA,mn_ou=NA,mx_ou=NA)
cnt=1
#ΔV / SA(T1) ≈ [PA(T2)·(SA(T2) − PA(T2))/P(T2) − PA(T1)·(SA(T1) − PA(T1))/P(T1)] / SA(T1)

for(i in 1:(nrow(vs)-1)){
  for(j in (i+1):nrow(vs)){
    hj=(vs$SA[j]-vs$PA[j])/vs$P[j]
    hi=(vs$SA[i]-vs$PA[i])/vs$P[i]
    Vj=vs$PA[j]*hj
    Vi=vs$PA[i]*hi
    delV=Vj-Vi
    
    tr$s1[cnt]=vs$VolumeScenario[i]
    tr$s2[cnt]=vs$VolumeScenario[j]
    tr$scen[cnt]=paste0(vs$VolumeScenario[i],vs$VolumeScenario[j])
    tr$delPA[cnt]=vs$PA[j]-vs$PA[i]
    tr$pctPA[cnt]=tr$delPA[cnt]/vs$PA[i]
    tr$delSA[cnt]=vs$SA[j]-vs$SA[i]
    tr$delP[cnt]=vs$P[j]-vs$P[i]
    tr$delV[cnt]=vs$V[j]-vs$V[i]
    tr$PlanarLE[cnt]=tr$delPA[cnt]/vs$P[i]
    tr$up[cnt]=(tr$delSA[cnt]-tr$delPA[cnt])/vs$PA[i]
    tr$delRug[cnt]=(vs$SA[j]/vs$PA[j])-(vs$SA[i]/vs$PA[i])
    tr$delH[cnt]=hj-hi
    tr$sap[cnt]=tr$delSA[cnt]/vs$P[i]
    tr$VolumeLE_proxy[cnt]=delV/vs$SA[i]
    tr$VolumeLE_check[cnt]=LE_vol(SA1 = vs$SA[i],PA1 = vs$PA[i],P1 = vs$P[i],
                               SA2 = vs$SA[j],PA2 = vs$PA[j],P2 = vs$P[j])
    tr$VolumeLE[cnt]=tr$delV[cnt]/vs$SA[i]
    tr$mn_ou[cnt]=mean(c(tr$PlanarLE[cnt],tr$up[cnt]))
    tr$mx_ou[cnt]=max(c(tr$PlanarLE[cnt],tr$up[cnt]))
    tr$mx_oh[cnt]=max(c(tr$PlanarLE[cnt],tr$delH[cnt]))
    cnt=cnt+1
  }
}
V=tr$VolumeLE
tr$scen


RMSE=function(x,y){
  E=x-y
  S=E^2
  M=mean(S)
  R=sqrt(M)
  return(R)
}

cor(V,tr$up)
cor(V,tr$mx_ou)
cor(V,tr$sap)
cor(V,tr$mn_ou)
cor(V,tr$delH)
cor(V,tr$mx_oh)
cor(V,tr$VolumeLE_proxy)
cor(V,tr$VolumeLE_check)
cor(V,tr$PlanarLE)

RMSE(V,tr$up)
RMSE(V,tr$mx_ou)
RMSE(V,tr$sap)
RMSE(V,tr$mn_ou)
RMSE(V,tr$delH)
RMSE(V,tr$mx_oh)
RMSE(V,tr$VolumeLE_proxy)
RMSE(V,tr$VolumeLE_check)
RMSE(V,tr$PlanarLE)

tr
library(tidyverse)
library(ggrepel)

tr %>% 
  #filter(!s1%in%c("K","M"),!s2%in%c("K","M")) %>% 
  ggplot(aes(x=VolumeLE,label=scen))+
  # geom_point(aes(y=mn_ou),color="gold",shape=1,size=4)+
  # geom_point(aes(y=mx_ou),color="gray",shape=1,size=4)+
  #geom_point(aes(y=VolumeLE_proxy),color="black",size=2)+
  geom_text(aes(y=VolumeLE_proxy,label=scen,color="VolumeLE"),size=3)+
  stat_smooth(aes(y=VolumeLE_proxy,color="VolumeLE"))+
  # geom_point(aes(y=VolumeLE_check),color="black",shape=1,size=4)+
  #  geom_point(aes(y=up),color="red")+
  #geom_point(aes(y=PlanarLE),color="blue")+
  geom_text(aes(y=PlanarLE+.001,label=scen,color="PlanarLE"),size=3)+
  stat_smooth(aes(y=PlanarLE+.001,color="PlanarLE"))+
  #  geom_point(aes(y=sap),color="darkgreen")+
  scale_color_manual(name="LE Metric",guide="legend",values=c("VolumeLE"="black","PlanarLE"="blue"))+
  scale_fill_viridis_c()+
  #geom_label_repel(aes(y=VolumeLE_proxy,fill=as.numeric(factor(pctPA))),alpha=.5)+
  geom_abline()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab("Actual Volumetric Linear Extension (cm)")+
  ylab("Planar and Volume Proxies Linear Extension (cm)")+
  scale_x_log10(limits=c(.01,20))+
  scale_y_log10(limits=c(.001,20))+
  theme_bw()+theme(legend.position = c(0.1,.9))+
  ggtitle("Volume Scenario Linear Extension Metrics - Planar & Volume Proxy VS Gold-Standard Volume")

tr %>% ggplot(aes(VolumeLE,mn_ou,label=scen))+geom_point()+geom_label_repel()+geom_abline()


plot(tr$up,tr$PlanarLE)
plot(tr$PlanarLE,tr$VolumeLE,col="red",pch=1,cex=2)

plot(tr$up,tr$delH,col="red",pch=1,cex=2)


points(tr$up,tr$VolumeLE,col="gold",pch=1,cex=2)
points(tr$delH,tr$VolumeLE,col="yellow",pch=1,cex=2)
points(tr$sap,tr$VolumeLE,col="blue",pch=16)
points(tr$mn_ou,tr$VolumeLE,col="darkgreen",pch=16)
abline(a=0,b=1)
cor(tr$sap,tr$VolumeLE)
cor(tr$mn_ou,tr$VolumeLE)

