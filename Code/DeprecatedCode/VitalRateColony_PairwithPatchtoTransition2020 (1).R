library(plyr)
library(lubridate)
library(ggplot2)

load(file = "/Users/c-rod/Documents/GitHub/Patch-To-Transition/ThesisSites_2021_09_CR.rdata")
VitalRate_Growth$TransitionType=factor(VitalRate_Growth$TransitionType,levels=c("RECR","MORT","GROWTH","SHRINK","FISSION","FUSION","FUSION_FISSION"))
#Aggregate Colonies
CID=substr(VitalRate_Growth$ColonyID,1,7)
CID2=paste0(substr(VitalRate_Growth$T0_PatchName,1,17),CID,substr(VitalRate_Growth$T0_PatchName,25,26))
RM_i=which(substr(CID2,1,4)%in%c("MORT","RECR"))
CID2[RM_i]=paste0(substr(VitalRate_Growth$T0_PatchName[RM_i],6,22),CID[RM_i],substr(VitalRate_Growth$T0_PatchName[RM_i],30,31))
VitalRate_Growth$ColonyID=CID2
VitalRate_Patch=subset(VitalRate_Growth,DataOrError=="DATA")  

head(VitalRate_Patch,5)

PatchSize=data.frame(Name=c(as.vector(VitalRate_Patch$T0_PatchName),as.vector(VitalRate_Patch$T1_PatchName)),
                     Size=c(VitalRate_Patch$StartingSize,VitalRate_Patch$EndingSize),
                     Perim=c(VitalRate_Patch$StartingPerimeter,VitalRate_Patch$EndingPerimeter),
                     Diam=c(VitalRate_Patch$StartingMaxDiam,VitalRate_Patch$EndingMaxDiam))

dim(PatchSize)
PatchSize=unique(PatchSize)
dim(PatchSize)
PatchSizeLU=PatchSize$Size
names(PatchSizeLU)=PatchSize$Name
PatchPerimLU=PatchSize$Perim
names(PatchPerimLU)=PatchSize$Name
PatchDiamLU=PatchSize$Diam
names(PatchDiamLU)=PatchSize$Name
VitalRate_Colony=ddply(VitalRate_Patch,.(Site,DataOrError,ColonyID,Spec_Code,Genus_Code,StartingDate,EndingDate,Interval_Years),summarize,
                       N_t0=length(unique(T0_PatchName)),
                       N_t1=length(unique(T1_PatchName)),
                       StartingSize=sum(PatchSizeLU[unique(as.vector(T0_PatchName))],na.rm=T),
                       EndingSize=sum(PatchSizeLU[unique(as.vector(T1_PatchName))],na.rm=T),
                       StartingPerim=sum(PatchPerimLU[unique(as.vector(T0_PatchName))],na.rm=T),
                       EndingPerim=sum(PatchPerimLU[unique(as.vector(T1_PatchName))],na.rm=T),
                       StartingSummedPatchDiam=sum(PatchDiamLU[unique(as.vector(T0_PatchName))],na.rm=T),
                       EndingSummedPatchDiam=sum(PatchDiamLU[unique(as.vector(T1_PatchName))],na.rm=T),
                       TransitionMagnitude=EndingSize-StartingSize,
                       TransitionType=paste0(unique(TransitionType),collapse="_"),
                       PercentChange=(EndingSize-StartingSize)/StartingSize,
                       Log2Ratio_Change=log2(EndingSize/StartingSize))
VitalRate_Colony$TransitionTypeSimple="GROWTH"
VitalRate_Colony$TransitionTypeSimple[VitalRate_Colony$StartingSize==0]="RECR"
VitalRate_Colony$TransitionTypeSimple[VitalRate_Colony$EndingSize==0]="MORT"
VitalRate_Colony$EndingPerim[VitalRate_Colony$TransitionTypeSimple=="MORT"]=0
VitalRate_Colony$EndingSummedPatchDiam[VitalRate_Colony$TransitionTypeSimple=="MORT"]=0
VitalRate_Colony$StartingPerim[VitalRate_Colony$TransitionTypeSimple=="RECR"]=0
VitalRate_Colony$StartingSummedPatchDiam[VitalRate_Colony$TransitionTypeSimple=="RECR"]=0
  
VitalRate_Colony$TransitionTypeSimple[VitalRate_Colony$TransitionTypeSimple=="GROWTH"&VitalRate_Colony$EndingSize<VitalRate_Colony$StartingSize]="SHRINK"
VitalRate_Colony$Fragmented=VitalRate_Colony$N_t0>1|VitalRate_Colony$N_t1>1

table(VitalRate_Colony$Fragmented,VitalRate_Colony$TransitionTypeSimple)
VitalRate_Colony[VitalRate_Colony$Fragmented==T&VitalRate_Colony$TransitionTypeSimple=="RECR",]


PlotGenusLU=c("Montipora sp.","Pocillopora sp.","Porites sp.");names(PlotGenusLU)=c("MOSP","POCS","POSP")
VitalRate_Colony$Genus=PlotGenusLU[as.vector(VitalRate_Colony$Genus_Code)]
VitalRate_Colony$Recruit=as.numeric(VitalRate_Colony$TransitionType=="RECR")
VitalRate_Colony$Mortality=as.numeric(VitalRate_Colony$TransitionType=="MORT")
VitalRate_Colony$Fragmentation=as.numeric(VitalRate_Colony$Fragmented)
VitalRate_Colony$ln_SS=log(VitalRate_Colony$StartingSize)
VitalRate_Colony$ln_ES=log(VitalRate_Colony$EndingSize)
VRCd=subset(VitalRate_Colony,DataOrError=="DATA")
VRCd$Genus=factor(VRCd$Genus)
VRCd$StartingYear=year(VRCd$StartingDate)
VRCd$EndingYear=year(VRCd$EndingDate)
VRCd$Interval=paste0(substr(VRCd$StartingYear,3,4),"-",substr(VRCd$EndingYear,3,4))
#VRCd$SiteInterval=paste0(substr(VRCd$Site,5,11),"\n",VRCd$Interval) #old line
VRCd$SiteInterval=paste0(substr(VRCd$Site,1,11),"\n",VRCd$Interval) #Changed SiteInterval to include island name in the site

VRCd$Island=factor(substr(VRCd$Site,1,3),levels=c("FFS","HAW","MAI","OAH","PHR","KUR"))
VRCd$PropMagnitude=VRCd$EndingSize/VRCd$StartingSize
VRCd$AnnualPropRate_E=(VRCd$PropMagnitude)^(1/VRCd$Interval_Years)
VRCd$TransitionRate_L=VRCd$TransitionMagnitude/VRCd$Interval_Years
VRCd$AnnualEndingSize_E=VRCd$StartingSize*VRCd$AnnualPropRate_E
VRCd$TriennialEndingSize_E=VRCd$StartingSize*(VRCd$AnnualPropRate_E^3)
VRCd$TriennialEndingSize_E[VRCd$TransitionTypeSimple=="MORT"]=0
VRCd$ln_AES=log(VRCd$AnnualEndingSize_E)
VRCd$ln_3ES=log(VRCd$TriennialEndingSize_E)
VRCd_gs=subset(VRCd,TransitionTypeSimple%in%c("GROWTH","SHRINK"))

PatchLevel=VitalRate_Patch
ColonyLevel=VRCd
save(list = c("PatchLevel","ColonyLevel"),file = "/Users/c-rod/Documents/GitHub/Patch-To-Transition/Patch_And_Colony_Data_20210917.rdata")


ggplot(VRCd,aes(x=StartingSize,y=AnnualPropRate_E,size=N_t0,shape=factor(Fragmentation)))+
  geom_point()+
  facet_grid(Genus~.)+
  geom_hline(yintercept = 1)+scale_x_sqrt(limits=c(0,2000))+
  scale_shape_manual(name="Fragmentation",values=c(21,4))+theme_bw()

ggplot(VRCd,aes(x=cut(N_t0,c(0,1.5,5,10,50,100)),y=AnnualPropRate_E))+
#  geom_point()+
  geom_boxplot()+scale_y_sqrt()+
  facet_grid(Genus~.)+
  geom_hline(yintercept = 1)+theme_bw()
  
VRCd$DeltaN=VRCd$N_t1-VRCd$N_t0

df=VRCd[,c("N_t0","N_t1")]
df=na.omit(df)
ggplot(VRCd,aes(x=N_t0,y=N_t1))+
  geom_jitter(colour="red",alpha=0.2,height=.25, width=.25)+
  theme_bw()
length(which(VRCd$N_t0<1))

ggplot(VRCd_gs,aes(x=SiteInterval,y=AnnualPropRate_E,fill=Island,color=Interval))+
  geom_boxplot()+
  facet_wrap(Genus_Code~.,scales="free")+
  theme_bw()+theme(axis.text=element_text(angle=90))

ggplot(subset(VRCd,Genus=="Porites sp."),aes(x=Island,y=AnnualPropRate_E,fill=Interval))+
  geom_boxplot(width=.2)

Pplot=ggplot(VitalRate_Patch,aes(log(StartingSize),log(EndingSize),color=TransitionType,shape=TransitionType))+geom_point(size=2)+geom_abline()+facet_grid(Genus_Code~.)+
  scale_shape_manual(values=c(3,4,1,0,6,2,11))+ggtitle("Patch Transition")+theme_bw()+theme(legend.position = "bottom")+xlab("Patch Size at T=0 (log(cm^2))")+ylab("Patch Size at T=1 (log(cm^2))")
Cplot=ggplot(VitalRate_Colony,aes(log(StartingSize),log(EndingSize),color=TransitionTypeSimple,shape=Fragmented))+geom_point(size=2)+geom_abline()+facet_grid(Genus_Code~.)+
  ggtitle("Colony Transition")+scale_shape_manual(values=c(1,4))+theme_bw()+theme(legend.position = "bottom")+xlab("Colony Size at T=0 (log(cm^2))")+ylab("Colony Size at T=1 (log(cm^2))")
library(gridExtra)
PC=grid.arrange(Pplot,Cplot,ncol=2)
sc=1
ggsave(plot=PC,filename = "M:/FixedSiteProducts/Vital Rates Data/Patch vs Colony.jpg",width=sc*16,height=sc*16)

CDplot=ggplot(VRCd,aes(StartingSize,EndingSize,fill=TransitionTypeSimple,color=TransitionTypeSimple,shape=Fragmented))+
  geom_point(size=2)+
  geom_abline()+
  facet_grid(.~Genus)+
  ggtitle("Colony Transition")+
  xlab("Colony Size at T=0 (cm^2)")+
  ylab("Colony Size at T=3 yrs (cm^2)")+
  scale_shape_manual(name="Fragmentation",values=c(21,4))+
  scale_fill_brewer(name="Transition Type",palette = "Set1")+
  scale_color_brewer(name="Transition Type",palette = "Set1")+
  scale_y_sqrt(limits=c(0,1500),breaks=round(c(0,1,2.5,5,10,15,20,25)^2*pi))+
  scale_x_sqrt(limits=c(0,1500),breaks=round(c(0,1,2.5,5,10,15,20,25)^2*pi))+
  coord_equal()+
  theme_bw()+theme(legend.position = "bottom")+
  guides(shape = guide_legend(override.aes = list(size = 5)),
         fill = guide_legend(override.aes = list(size = 5)),
         color = guide_legend(override.aes = list(size = 5)))
CDplot

sc=2.5
ggsave(plot = CDplot,
       filename = paste0("C:/Users/Thomas.Oliver/WORK/Projects/Vital Rates/ErrorRates/Colony_T0T1Size.jpg"),
       width=sc*6,height=sc*3)

  
# Plot Recruit Rate ####
library(ggpubr)
library(viridis)

RecDen=ddply(VRCd,.(Site,EndingDate,Spec_Code,Genus),summarize,
               Nrecruit=sum(Recruit))
names(RecDen)[2]="Date"
Patches$Genus=PlotGenusLU[as.vector(Patches$Genus_Code)]
SA=ddply(Patches,.(Site,Date,Spec_Code),summarize,Ncirc=max(Circrat),Area=Ncirc*2.5)
SA$Date=mdy(SA$Date)
RecDen=merge(RecDen,SA)
RecDen$Nrec_area=RecDen$Nrecruit/RecDen$Area
RecDenG=ddply(RecDen,.(Site,Date,Genus),summarize,Nrec_area=sum(Nrec_area))

RecDenPlot=ggplot(RecDenG,aes(x=Genus,y=Nrec_area,fill=Genus))+
  geom_violin(width=0.9,scale="width",alpha=0.33)+
  geom_boxplot(width=0.33,color="black",size=1)+
  scale_fill_viridis(discrete = T,end=.85)+
  xlab("Genus")+
  ylab("Recruits per Area (N/m^2)")+
  ggtitle("Genus-Specific Recruit Density")+
  theme_bw()+
  theme(legend.position = "none")
scale=1
ggsave(filename = paste0("C:/Users/Thomas.Oliver/WORK/Projects/Vital Rates/ErrorRates/Colony_RecruitDensity.jpg"),plot=RecDenPlot,width=6*scale*1.1618,height=scale*6)

# Plot Size Specific Mortality Rate ####
modMort=glm(formula=Mortality~StartingSize:Genus,data=VRCd,family=binomial)
summary(modMort)

#Prediction for each Site, Genus, Size Combo
SizeRange=exp(seq(0,log(1500),length.out=100))
uGenera=unique(VRCd$Genus)
nSizeReps=length(uGenera)
nGenReps=length(SizeRange)
pred_VRCd=data.frame(Genus=sort(rep(uGenera,nGenReps)),
                          StartingSize=rep(SizeRange,nSizeReps))
#subset predictions within size range of data...
uG=unique(pred_VRCd$Genus)
predsub=NULL
SR=ddply(VRCd,.(Genus),summarize,
         Size_min=quantile(StartingSize,0,na.rm=T),
         Size_q01=quantile(StartingSize,.01,na.rm=T),
         Size_q99=quantile(StartingSize,.99,na.rm=T),
         Size_max=quantile(StartingSize,1,na.rm=T))


for(j in 1:length(uG)){
  this_SG=subset(pred_VRCd,
                 Genus==uG[j]&
                   #StartingSize>subset(SR,Genus==uG[j])$Size_q01&
                   StartingSize<subset(SR,Genus==uG[j])$Size_q99)
  predsub=rbind(predsub,this_SG)
}


pModMort=predict(modMort,se=T,
                 newdata=predsub,
                 type="response")
predsub$Mortality=pModMort$fit
predsub$Mort_min=pModMort$fit-pModMort$se.fit
predsub$Mort_max=pModMort$fit+pModMort$se.fit

SSM=ggplot()+
  geom_vline(xintercept=pi*(2.5^2),color="gray50")+
  geom_hline(yintercept=c(0,1),color="gray50")+
  geom_jitter(data=VRCd,aes(x=StartingSize,y=Mortality,color=Genus),height=.05)+
  geom_ribbon(data=predsub,aes(x=StartingSize,ymin=Mort_min,ymax=Mort_max,fill=Genus),alpha=.1)+
  geom_line(data=predsub,aes(x=StartingSize,y=Mortality,color=Genus),size=1.25,lty=1)+
  ylab("Mortality (Prob.)")+
  xlab("Starting Colony Area (cm2)")+
  ggtitle("Size-Dependent Mortality")+
  scale_color_viridis(discrete = T,end=.85)+
  scale_fill_viridis(discrete = T,end=.85)+
  facet_grid(.~Genus)+
  scale_x_log10(limits=c(1,1500),breaks=round(c(1,2.5,5,10,15,20,25)^2*pi))+
  theme_bw()+
  theme(legend.position = "bottom")+
  guides(shape = guide_legend(override.aes = list(size = 5)),
         fill = guide_legend(override.aes = list(size = 5)),
         color = guide_legend(override.aes = list(size = 5)))
scale=2.5
SSM
ggsave(filename = paste0("C:/Users/Thomas.Oliver/WORK/Projects/Vital Rates/ErrorRates/Colony_Mortality.jpg"),
       plot=SSM,width=6*scale,height=scale*3)




# Plot Size Specific Fragmentation Rate ####
modFrag=glm(formula=Fragmentation~StartingSize:Genus,data=VRCd,family=binomial)
summary(modFrag)


pModFrag=predict(modFrag,se=T,
                 newdata=predsub,
                 type="response")
predsub$Fragmentation=pModFrag$fit
predsub$Frag_min=pModFrag$fit-pModFrag$se.fit
predsub$Frag_max=pModFrag$fit+pModFrag$se.fit

SSF=ggplot()+
  geom_vline(xintercept=pi*(2.5^2),color="gray50")+
  geom_hline(yintercept=c(0,1),color="gray50")+
  geom_jitter(data=VRCd,aes(x=StartingSize,y=Fragmentation,color=Genus),height=.05)+
  geom_ribbon(data=predsub,aes(x=StartingSize,ymin=Frag_min,ymax=Frag_max,fill=Genus),alpha=.1)+
  geom_line(data=predsub,aes(x=StartingSize,y=Fragmentation,color=Genus),size=1.25,lty=1)+
  ylab("Fragmentation (Prob.)")+
  xlab("Starting Colony Area (cm2)")+
  ggtitle("Size-Dependent Fragmentation")+
  scale_fill_viridis(discrete = T,end=.85)+
  scale_color_viridis(discrete = T,end=.85)+
  facet_grid(.~Genus)+
  scale_x_log10(limits=c(1,1500),breaks=round(c(1,2.5,5,10,15,20,25)^2*pi))+
  theme_bw()+
  theme(legend.position = "bottom")+
  guides(shape = guide_legend(override.aes = list(size = 5)),
         fill = guide_legend(override.aes = list(size = 5)),
         color = guide_legend(override.aes = list(size = 5)))

SSF
scale=2.5
ggsave(filename = paste0("C:/Users/Thomas.Oliver/WORK/Projects/Vital Rates/ErrorRates/Colony_Fragmentation.jpg"),
       plot=SSF,width=6*scale,height=scale*3)


#Size Specific Growth Rate
library(mgcv)
VRCd_gsUF=subset(VRCd_gs,Fragmented==F)
q80=quantile(VRCd_gsUF$TransitionRate,.8)
VRCd_gsUF80=subset(VRCd_gsUF,TransitionRate>=q80)
#modGrow=glm(formula=ln_3ES~ln_SS:Genus,data=VRCd_gs,family=gaussian)
head(VRD_d)
modGrow_0=gam(formula=AnnualEndingSize_E~1,data=VRCd_gs)
modGrow_1.1=gam(formula=AnnualEndingSize_E~Genus,data=VRCd_gs)
modGrow_1.2=gam(formula=AnnualEndingSize_E~s(StartingSize,by=Genus,k=3),data=VRCd_gs)
modGrow_1.3=gam(formula=AnnualEndingSize_E~Genus+s(StartingSize,by=Genus,k=3),data=VRCd_gs)
modGrow_2=gam(formula=AnnualEndingSize_E~s(StartingSize,by=Genus,k=3)+s(N_t0,by=Genus,k=3),data=VRCd_gs)
modGrow_3.10=gam(formula=AnnualEndingSize_E~s(StartingSize,by=Genus,k=3)+s(N_t0,by=Genus,k=3)+Interval,data=VRCd_gs)
modGrow_3.11=gam(formula=AnnualEndingSize_E~s(StartingSize,by=Genus,k=3)+s(N_t0,by=Genus,k=3)+Site,data=VRCd_gs)
modGrow_3.2=gam(formula=AnnualEndingSize_E~s(StartingSize,by=Genus,k=3)+s(N_t0,by=Genus,k=3)+Site+Interval,data=VRCd_gs)
modGrow_3.3=gam(formula=AnnualEndingSize_E~s(StartingSize,by=Genus,k=3)+s(N_t0,by=Genus,k=3)+SiteInterval,data=VRCd_gs)

mGrowth=glm(formula=log(AnnualEndingSize_E)~Genus:log(StartingSize)+N_t0+Site+Interval,data=VRCd_gs)
mGrowth2=glm(formula=log(AnnualEndingSize_E)~Genus:log(StartingSize)+N_t0+SiteInterval,data=VRCd_gs)
mMort=glm(formula=Mortality~StartingSize:Genus+Site+Interval,data=subset(VRCd,TransitionTypeSimple!="RECRUIT"),family=binomial)
sGro=summary(mGrowth)
sGro2=summary(mGrowth2)
rownames(sGro$coefficients)[3:19]
rownames(sGro2$coefficients)[3:19]

library(sjPlot)
plot_model(mGrowth,terms = rownames(sGro$coefficients)[3:19])+geom_vline(xintercept = 0,color="black")+ggtitle("Colony Growth: Effect of Site and Interval")
plot_model(mGrowth2,terms = rownames(sGro2$coefficients)[3:19])+geom_vline(xintercept = 0,color="black")+ggtitle("Colony Growth: Effect of Site and Interval")
plot_model(mMort)


BIC(modGrow_0,modGrow_1.1)
BIC(modGrow_0,modGrow_1.2)
BIC(modGrow_1.1,modGrow_1.2)
BIC(modGrow_1.1,modGrow_1.3)
BIC(modGrow_1.2,modGrow_1.3)
BIC(modGrow_1.2,modGrow_2)
BIC(modGrow_2,modGrow_3.10)
BIC(modGrow_3.10,modGrow_3.11)
BIC(modGrow_3.10,modGrow_3.2)
BIC(modGrow_3.10,modGrow_3.3)
AIC(modGrow_3.10,modGrow_3.2)
AIC(modGrow_3.2,modGrow_3.3)

mod_Grow=modGrow_3.2
summary(modGrow)
plot(modGrow,pages=1,residuals = T)

ggplot(VRCd_gs,aes(StartingSize,AnnualEndingSize_E,color=Site))+
  facet_grid("Genus")+
  geom_point()+geom_smooth(span=2)

predsub$ln_SS=log(predsub$StartingSize)
pModGrow=predict(modGrow,se=T,
                 newdata=predsub,
                 type="response")
predsub$Growth=pModGrow$fit
predsub$Growth_min=pModGrow$fit-pModGrow$se.fit
predsub$Growth_max=pModGrow$fit+pModGrow$se.fit
#predsub$eGrowth=exp(predsub$Growth)
#predsub$eGrowth_min=exp(predsub$Growth_min)
#predsub$eGrowth_max=exp(predsub$Growth_max)

GbyS=ggplot()+
  geom_jitter(data=VRCd_gsUF,aes(x=StartingSize,y=AnnualEndingSize,color=Genus),size=1,width=.1,height=.1,alpha=.5)+
  geom_ribbon(data=predsub,aes(x=StartingSize,ymin=Growth_min,ymax=Growth_max,fill=Genus),alpha=.75)+
  geom_line(data=predsub,aes(x=StartingSize,y=Growth,color=Genus),size=1)+
  ylab("Colony Area + 1 yr (cm2)")+
  xlab("Colony Area (cm2)")+
  ggtitle("Size Dependent Growth")+
  #geom_hline(yintercept=0,color="black")+
  geom_abline()+
  geom_vline(xintercept=pi*(2.5^2),color="gray50")+
  geom_hline(yintercept=pi*(2.5^2),color="gray50")+
  facet_grid(.~Genus)+
  scale_fill_viridis(discrete = T,end=.85)+
  scale_color_viridis(discrete = T,end=.85)+
  scale_y_log10(limits=c(1,1500),breaks=round(c(1,2.5,5,10,15,20,25)^2*pi))+
  scale_x_log10(limits=c(1,1500),breaks=round(c(1,2.5,5,10,15,20,25)^2*pi))+
  coord_equal()+
  theme_bw()+
  theme(legend.position = "bottom")+
  guides(shape = guide_legend(override.aes = list(size = 5)),
         fill = guide_legend(override.aes = list(size = 5)),
         color = guide_legend(override.aes = list(size = 5)))
GbyS
scale=2.5
ggsave(filename = paste0("C:/Users/Thomas.Oliver/WORK/Projects/Vital Rates/ErrorRates/Colony_Growth.jpg"),
       plot=GbyS,width=6*scale,height=scale*3)

VRCd$modAnnualEndingSize=predict(modGrow,newdata = data.frame(VRCd[,c("Genus","StartingSize")]))
VRCd$modPartialMortality=VRCd$modAnnualEndingSize-VRCd$AnnualEndingSize
VRCd$modPartialMortality[VRCd$modPartialMortality<0]=0
VRCd$modPartialMortality[VRCd$TransitionTypeSimple=="MORT"]=VRCd$StartingSize[VRCd$TransitionTypeSimple=="MORT"]
VRCd$modPctPartialMortality=VRCd$modPartialMortality/VRCd$StartingSize
ggplot(VRCd,aes(StartingSize,modPctPartialMortality))+
  geom_point()+
  facet_wrap("Genus",scales = "free")


AdSize=2.5^2*pi

PropMort=ddply(subset(VRCd,StartingSize>AdSize),.(Site,Spec_Code,Genus,StartingDate,EndingDate),summarize,
               Ntrans=length(StartingSize),
               pM=sum(Mortality)/length(Mortality),
               meanG=mean(100*PercentChange,na.rm=T),
               seG=sd(100*PercentChange,na.rm=T)/sqrt(length(PercentChange)))

RecDen=ddply(VRCd,.(Site,Spec_Code,Genus,EndingDate),summarize,
             Nrecruit=sum(Recruit))
names(RecDen)[4]="Date"
Patches$Genus=PlotGenusLU[as.vector(Patches$Genus_Code)]
SA=ddply(Patches,.(Site,Date,Spec_Code,Genus),summarize,Ncirc=max(Circrat),Area=Ncirc*2.5)
SA$Date=mdy(SA$Date)
RecDen=merge(RecDen,SA)
RecDen$Nrec_area=RecDen$Nrecruit/RecDen$Area
names(RecDen)[4]="EndingDate"
dim(PropMort)
PropMort=merge(PropMort,RecDen)
dim(PropMort)

PropMort$Island=substr(PropMort$Site,1,3)
PropMort$SiteNum=substr(PropMort$Site,9,11)
PropMort$YrInt=paste0((year(PropMort$StartingDate)-2000),"-",(year(PropMort$EndingDate)-2000))
PropMort$StTP=paste0(PropMort$SiteNum,"\n",PropMort$YrInt)
PropMort$Island=factor(PropMort$Island,levels=c("HAW","MAI","OAH","PHR","KUR"))
PropMort$Genus=factor(PropMort$Genus,levels=c("Pocillopora sp.","Porites sp.","Montipora sp."))
PropMort$Spec_Code=factor(PropMort$Spec_Code,levels=c("PMEA","PGRA","POCS","PLOB","PLUT","PLIG","POSP","MCAP","MPAT","MFLA","MOSP"))

library(ggrepel)
ggplot(PropMort,aes(x=pM,y=meanG,ymin=meanG-seG,ymax=meanG+seG,shape=Spec_Code,color=Genus))+
  geom_point(aes(size=Nrec_area))+
  geom_errorbar()+
  geom_text_repel(aes(label=StTP),size=3)+
  facet_wrap("Island")+
  geom_hline(yintercept=0)+
  scale_shape_manual(values=0:10)+
  scale_size(name="Recruit Density (N/m2)")+
  theme_bw()+xlab("Proportion Whole Colony Mortality")+ylab("Mean Percentage Growth")


################################################################################################################################
#
JuvA=2.5^2*pi
VitalRate_Patch$AdJuv="AD"
VitalRate_Patch$AdJuv[VitalRate_Patch$EndingSize<JuvA]="JUV"
VitalRate_Patch$AdJuv[VitalRate_Patch$EndingSize==0]="MORT"
VitalRate_Patch$TrueRec="OTHER"
VitalRate_Patch$TrueRec[VitalRate_Patch$TransitionType=="RECR"]="RECR"

rectab=table(VitalRate_Patch$AdJuv,VitalRate_Patch$TrueRec,VitalRate_Patch$Genus_Code)
TruePos=100*rectab[2,2,]/(rectab[2,1,]+rectab[2,2,])
FalsePos=100*rectab[2,1,]/(rectab[2,1,]+rectab[2,2,])

rectab
TruePos
FalsePos

JuvOnly=subset(VitalRate_Patch,AdJuv=="JUV")
table(JuvOnly$TransitionType,JuvOnly$Genus_Code)
