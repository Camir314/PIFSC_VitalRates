#First up, pool "REPLACEMENT" SIGs and run TRUE/FALSE g/S models by genus
#rm(list=ls())
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
library(gridExtra)
library(boot)
library(ggrepel)
library(ggbreak)
# Loading Data -------------------------------------------------------
loadnames1=load("./Data/ColonyTransitions/Script_Step3_DataPackage.rdata");loadnames1
loadnames2=load("./Data/ModelData/HA_MA_AS_SIG_IMP_OUTPUTS.rdata");loadnames2
loadnames3=load("./Data/ColonyTransitions/HA_MA_AS_ColonyTransition_withAnnSurv.rdata");loadnames3

#First Pull Together Models from Current and Former Work
# Global mesh variables 
MatrixVal=list(
  min.size = -0.5,   # 2*sqrt((10^-2.5)/pi)  # approx 6 cm  diameter coral
  max.size = 1.1*max(c(ColTrans$log10_SS, ColTrans$log10_ES), na.rm = T),
  n = 50, #mesh size / number of cells in the discretized kernel
  rec.size = -0.1)#,  # 2*sqrt((10^-0.1)/pi)  # approx 1 cm  diameter coral #2*sqrt((10^1.29303)/pi)  # approx 5 cm  diameter coral
MatrixVal$bin_size = MatrixVal$min.size + c(0:MatrixVal$n) * (MatrixVal$max.size - MatrixVal$min.size)/MatrixVal$n #boundary points (the edges of the cells defining the kernel 
MatrixVal$y = 0.5 * (MatrixVal$bin_size[1:MatrixVal$n]+MatrixVal$bin_size[2:(MatrixVal$n+1)]) #mesh points (midpoints of cells)
MatrixVal$I = MatrixVal$y >= MatrixVal$rec.size
MatrixVal$delta_size = MatrixVal$y[2] - MatrixVal$y[1] #width of cells (h)

s2sc=Site2Sec %>% st_coordinates() %>% as.data.frame() %>% rename(Longitude=X,Latitue=Y) 
Site2SecLL=Site2Sec %>% st_drop_geometry() %>% cbind(s2sc)%>% select(!c(Date,SEC_NAME,Region)) %>% distinct()
MODs=MODs %>% left_join(Site2SecLL,by=c("Site"))

#ALE
ColTrans$AnnualLinearExtension.cm[is.infinite(ColTrans$AnnualLinearExtension.cm)]=NA
Mean_SIG_ALE=ColTrans %>% 
  #filter(!is.na(SIG_REPLACEMENT)) %>% 
  group_by(SIG) %>% 
  summarise(ALE=mean(AnnualLinearExtension.cm,na.rm=T))
MODs=MODs %>% left_join(Mean_SIG_ALE,by="SIG")
MODs$SIG=factor(MODs$SIG,levels=MODs$SIG[order(MODs$Lambda,decreasing = T)])


#REPLACEMENT TRUE Model
GrowthR2_Standard <- 0.5 
MinN_Standard <- 6
# ColTrans=ColTrans %>% left_join(MODs[,c("SIG","REPLACEMENT")],by="SIG") %>% 
#   rename("SIG_REPLACEMENT"="REPLACEMENT")

PlanarGrowthReplacementTib=ColTrans %>% 
  filter(TransitionType %in% c("GROW"),!is.na(SIG_REPLACEMENT)) %>% 
  group_by(GENUS_CODE,SIG_REPLACEMENT) %>% 
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
head(PlanarGrowthReplacementTib)
PlanarGrowthReplacementDF = as.data.frame(PlanarGrowthReplacementTib)

### JUST REPLACEMENT by Genus
ReplacementAnnSurvTib=ColTrans %>% 
  filter(TransitionType %in% c("GROW","MORT")&!is.na(log10_SS)&!is.na(SIG_REPLACEMENT)) %>% 
  filter(s.BadModel==FALSE) %>% 
  group_by(GENUS_CODE,SIG_REPLACEMENT) %>% 
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
    s.BadModel=ifelse((is.na(s.AUC)|(s.AUC<0.50)|(s.AUC==1)),TRUE,FALSE),
  )
ReplacementAnnSurvDataFrame = as.data.frame(ReplacementAnnSurvTib)
(ReplacementAnnSurvDataFrame)


#Recruits
#rec values by region and genus code

#Replacement models (genus and region)
RecSecRepSummary=MODs %>%
  group_by(GENUS_CODE,REPLACEMENT) %>% 
  summarise(
    SecRecVal=median(SecRecVal)
  )
MODs %>% ggplot(aes(x=SecRecVal,fill=SecRecValSource))+
  geom_histogram(binwidth=0.05)+
  geom_vline(aes(xintercept=SecRecVal),data=RecSecRepSummary,fill="red")+
  scale_x_log10()

#add SfM sector-level stock recruitment data
combo_Replacement <- left_join(PlanarGrowthReplacementDF,ReplacementAnnSurvDataFrame)#,SAGrowthRegionalDF)
combo_Replacement <- (combo_Replacement)%>% filter(pg.BadModel==FALSE&s.BadModel==FALSE)

ModelParams_Replacement = combo_Replacement %>% 
  left_join(RecSecRepSummary, by= c("GENUS_CODE","SIG_REPLACEMENT"="REPLACEMENT") ) %>% 
  arrange(SIG_REPLACEMENT,GENUS_CODE)

Name="HA_MA_AS"
save(ModelParams_Replacement, file = sprintf("./Data/ModelData/%s_ReplacementModFits_noBADs.rdata",Name))

#compare SIGs with replacement
uG=unique(ModelParams_Replacement$GENUS_CODE)
delGen=data.frame(GENUS_CODE=uG)
y=MatrixVal$y
for(i in 1:length(uG)){
  ThisModT=ModelParams_Replacement %>% filter(GENUS_CODE==uG[i],SIG_REPLACEMENT==T)
  ThisModF=ModelParams_Replacement %>% filter(GENUS_CODE==uG[i],SIG_REPLACEMENT==F)
  GyT=ThisModT$pg.int+ThisModT$pg.slp*y
  GyF=ThisModF$pg.int+ThisModF$pg.slp*y
  SyT=inv.logit(ThisModT$s.int+ThisModT$s.slp*y)
  SyF=inv.logit(ThisModF$s.int+ThisModF$s.slp*y)
  # par(mfrow=c(2,2))
  # plot(y,GyT,col="darkgreen",type="l")
  # points(y,GyF,col="red",type="l")
  # plot(y,SyT,col="darkgreen",type="l")
  # points(y,SyF,col="red",type="l")
  # plot(ThisModT$SecRecVal,ThisModF$SecRecVal,col="darkgreen",type="p")
  # abline(0,1)
  #   delGen$dG[i]=mean((GyF-GyT)/GyT)
  delGen$dS[i]=mean((SyF-SyT)/SyT)
  delGen$dR[i]=mean((ThisModF$SecRecVal-ThisModT$SecRecVal)/ThisModT$SecRecVal)
}
delGen

y=MatrixVal$y
for(i in 1:nrow(MODs)){
  RepMatch=ModelParams_Replacement %>%
    filter(GENUS_CODE==MODs[i,"GENUS_CODE"]&SIG_REPLACEMENT==TRUE)
  Gy.=MODs$pg.int[i]+MODs$pg.slp[i]*y
  GyT=RepMatch$pg.int+RepMatch$pg.slp*y
  Sy.=inv.logit(MODs$s.int[i]+MODs$s.slp[i]*y)
  SyT=inv.logit(ThisModT$s.int+ThisModT$s.slp*y)
  MODs$delG[i]=mean((Gy.-GyT))
  MODs$delS[i]=mean((Sy.-SyT))
  MODs$delR[i]=mean((MODs$SecRecVal[i]-RepMatch$SecRecVal))
}

x1end=1.5
x2start=3
p=ggplot(MODs,aes(ALE,delS,size=SecRecVal,fill=Lambda,shape=REPLACEMENT))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  #facet_grid(GENUS_CODE~.)+
  scale_fill_viridis_c()+
  scale_shape_manual(values = c(21,22))+
  theme_bw()#+
#xlim(c(-.75,3.25))+
#scale_x_break(c(1.25,3))

MODpca=prcomp(x=scale(MODs[,c("SecRecVal","ALE","delS")]))
plot(MODpca)
Load=as.data.frame(MODpca$rotation);Load$Name=c("Recruitment","Growth","Survivorship")

MODs_pca=cbind(MODs,as.data.frame(MODpca$x))
MODs_pca$REPLACE_B=as.character(MODs_pca$REPLACEMENT)
MODs_pca$REPLACE_B[MODs_pca$REPLACE_B=="TRUE"]="ABOVE"
MODs_pca$REPLACE_B[MODs_pca$REPLACE_B=="FALSE"]="BELOW"
MODs_pca$REPLACE_B[MODs_pca$Lambda<1&MODs_pca$Lambda>=0.95]="BORDERLINE"
MODs_pca$REPLACE_B=factor(MODs_pca$REPLACE_B,levels=c("BELOW","BORDERLINE","ABOVE"))

ProjMat=matrix(c(1,1,-1,1,0,-1),ncol=2,byrow=T);
rownames(ProjMat)=c("G","S","R")
colnames(ProjMat)=c("ProjX","ProjY")
ProjDF=as.data.frame(ProjMat);ProjDF$Name=c("Growth","Survival","Recruitment")
MODS_scale=scale((MODs[,c("ALE","delS","SecRecVal")]),center = T,scale = T)
MODs_proj=MODS_scale%*%ProjMat
MODs_proj=as.data.frame(MODs_proj);names(MODs_proj)=c("ProjX","ProjY")
MODs_pca=cbind(MODs_pca,MODs_proj)

REPREF=MODs %>% filter(REPLACEMENT==T) %>%  group_by(GENUS_CODE) %>% summarize(
  ALE=mean(ALE),delS=mean(delS),SecRecVal=mean(SecRecVal)
)
REPREFscp=data.frame(GENUS_CODE=REPREF$GENUS_CODE,scale(REPREF[,2:4],center=attr(MODS_scale,"scaled:center"),scale=attr(MODS_scale,"scaled:scale"))%*%ProjMat)

Load_scale=2
MODs_pca %>% ggplot()+
  geom_point(aes(x=PC1,y=PC2,fill=Lambda,shape=REPLACE_B),size=3)+
  geom_segment(aes(xend=Load_scale*PC1,x=0,yend=Load_scale*PC2,y=0),
               arrow = arrow(length = unit(0.5,"cm")),
               data=Load)+
  geom_text_repel(aes(x=Load_scale*PC1,y=Load_scale*PC2,label=Name),data=Load)+
  scale_fill_viridis_c(breaks=c(.25,.5,.75,1),name="Pop. Growth")+
  scale_shape_manual(values = c(21,22,24),name="Pop. Replacement")+
  facet_wrap(.~GENUS_CODE)+
  theme_bw()


axes_range=1.1*max(abs(c(MODs_pca$ProjX,MODs_pca$ProjY)))
PROJplot=MODs_pca %>% ggplot()+
  geom_point(aes(x=ProjX,y=ProjY,fill=Lambda,shape=REPLACE_B,size=REPLACE_B))+
  geom_point(aes(x=ProjX,y=ProjY),fill="blue",shape=24,size=5,data=REPREFscp)+
  geom_segment(aes(xend=Load_scale*ProjX,x=0,yend=Load_scale*ProjY,y=0),
               arrow = arrow(length = unit(0.5,"cm")),
               data=ProjDF)+
  geom_text_repel(aes(x=1.2*Load_scale*ProjX,y=1.2*Load_scale*ProjY,label=Name),data=ProjDF)+
  scale_fill_viridis_c(breaks=c(.25,.5,.75,1),name="Pop. Growth")+
  scale_shape_manual(values = c(21,22,24),name="Pop. Replacement")+
  facet_wrap(.~GENUS_CODE,ncol=4,scales="free")+#xlim(c(-axes_range,axes_range))+ylim(c(-axes_range,axes_range))+
  theme_bw()+#coord_equal()+
  #  theme(legend.position = "bottom",legend.direction = "horizontal",legend.box = "vertical")
  theme(legend.position = "right",legend.direction = "vertical",legend.box = "vertical")+
  xlab("Projected Axis X")+ylab("Projected Axis Y")
PROJplot
PROJplot+xlim(c(-4,4))+ylim(c(-4,4))


ggplot(MODs_pca,aes(x=ALE,y=delS,shape=REPLACE_B,size=sqrt(SecRecVal),fill=Lambda))+
  geom_point()+scale_fill_viridis_c(breaks=c(.25,.5,.75,1),name="Pop. Growth")+
  scale_shape_manual(values = c(21,22,24),name="Pop. Replacement")+
  facet_wrap(REGION~GENUS_CODE)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)

ggplot(MODs_pca,aes(x=SecRecVal,y=delS,shape=REPLACE_B,size=ALE,fill=Lambda))+
  geom_point()+scale_fill_viridis_c(breaks=c(.25,.5,.75,1),name="Pop. Growth")+
  scale_shape_manual(values = c(21,22,24),name="Pop. Replacement")+
  facet_grid(REGION~GENUS_CODE)+geom_hline(yintercept = 0)+
  theme_bw()+scale_x_sqrt()


lenNONA=function(x){return(length(x[!is.na(x)]))}
se=function(x){return(sd(x,na.rm=T)/sqrt(lenNONA(x)))}
Grb=MODs_pca %>% group_by(GENUS_CODE,REPLACE_B) %>% summarize(ALEmn=mean(ALE),ALEse=se(ALE)) %>% mutate(YMIN=0,YMAX=18)
Gp=MODs_pca %>% ggplot()+
  # geom_rect(aes(xmin=ALEmn-ALEse,xmax=ALEmn+ALEse,
  #               ymin=YMIN,ymax=YMAX,
  #               fill=REPLACE_B),alpha=0.25,data=Grb,inherit.aes = F)+
  geom_vline(aes(xintercept=ALEmn,color=REPLACE_B),data=Grb,lty=2,size=1)+
  geom_histogram(aes(x=ALE,fill=REPLACE_B),binwidth=0.33)+
  facet_wrap(.~GENUS_CODE,scale="free_y",ncol=4)+
  geom_vline(xintercept = 0)+
  xlab("Growth - SIG Mean Annual Linear Extension (cm)")+theme_bw()+
  scale_color_discrete(guide="none")+
  scale_fill_discrete(name="Pop. Replacement:")+
  theme(legend.position = "bottom")+
  ggtitle("Vital Rates at 93 Sites-Interval-Genus Combos")
Srb=MODs_pca %>% group_by(GENUS_CODE,REPLACE_B) %>% summarize(delSmn=mean(delS),delSse=se(delS)) %>% mutate(YMIN=0,YMAX=18)
Sp=MODs_pca %>% ggplot(aes(x=delS,fill=REPLACE_B))+
  # geom_rect(aes(xmin=delSmn-delSse,xmax=delSmn+delSse,
  #               ymin=YMIN,ymax=YMAX,
  #               fill=REPLACE_B),alpha=0.25,data=Srb,inherit.aes = F)+
  geom_vline(aes(xintercept=delSmn,color=REPLACE_B),data=Srb,lty=2,size=1)+
  geom_histogram(binwidth=0.1)+geom_vline(xintercept = 0)+
  facet_wrap(.~GENUS_CODE,scale="free_y",ncol=4)+
  xlab("Survival\n(% Difference from Replacement Reference)")+theme_bw()+
  scale_fill_discrete(name="Pop. Replacement:")+
  scale_color_discrete(guide="none")+
  theme(legend.position = "bottom")
Rrb=MODs_pca %>% group_by(GENUS_CODE,REPLACE_B) %>% summarize(Recmn=mean(SecRecVal),Recse=se(SecRecVal)) %>% mutate(YMIN=0,YMAX=18)
Rp=MODs_pca %>% ggplot(aes(x=SecRecVal,fill=REPLACE_B))+
  # geom_rect(aes(xmin=Recmn-Recse,xmax=Recmn+Recse,
  #               ymin=YMIN,ymax=YMAX,
  #               fill=REPLACE_B),alpha=0.25,data=Rrb,inherit.aes = F)+
  geom_vline(aes(xintercept=Recmn,color=REPLACE_B),data=Rrb,lty=2,size=1)+
  geom_histogram(binwidth=0.01,aes(alpha=SecRecValSource))+geom_vline(xintercept = 0)+
  facet_wrap(.~GENUS_CODE,scale="free_y",ncol=4)+
  xlab("Recruitment\n(Sector-Level Recruits Per cm2 Adult Area)")+theme_bw()+
  scale_x_sqrt()+scale_alpha_manual(name="Recruit Estimate",labels=c("Sector","Region Ref."),values=c(1,.5))+
  scale_fill_discrete(name="Pop. Replacement:")+
  scale_color_discrete(guide="none")+
  theme(legend.position = "bottom")
VRHist=Gp/Sp/Rp+plot_layout(guides = "collect")&theme(legend.position = "bottom")
PROJplot/VRHist

MODs_pca=MODs_pca %>% arrange(REGION,Island,Site,Interval,GENUS_CODE)
SiteRank=MODs_pca %>% group_by(Site) %>% summarize(Lp=prod(Lambda)) %>% arrange(desc(Lp))
MODs_pca$Site=factor(MODs_pca$Site,levels=SiteRank$Site)
GC="MOSP"
Rc=MODs_pca %>% filter(GENUS_CODE==GC) %>% ggplot(aes(x=Site,y=SecRecVal,fill=REPLACE_B))+
  geom_col(position="dodge")+theme_bw()+facet_grid(.~REGION,scales = "free",drop=T)+
  geom_vline(xintercept = 0:99+.5,col="black",lty=2)+scale_y_sqrt()
Gc=MODs_pca %>% filter(GENUS_CODE==GC) %>% ggplot(aes(x=Site,y=ALE,fill=REPLACE_B))+
  geom_col(position="dodge")+theme_bw()+facet_grid(.~REGION,scales = "free",drop=T)+
  geom_vline(xintercept = 0:99+.5,col="black",lty=2)+geom_hline(yintercept = 0)
Sc=MODs_pca %>% filter(GENUS_CODE==GC) %>% ggplot(aes(x=Site,y=delS,fill=REPLACE_B))+
  geom_col(position="dodge")+theme_bw()+facet_grid(.~REGION,scales = "free",drop=T)+
  geom_vline(xintercept = 0:99+.5,col="black",lty=2)+geom_hline(yintercept = 0)
Gc/Sc/Rc

uSIG=unique(MODs$SIG)
#ROS-010_1807-2307_POSP
for (i in 1:length(uSIG)){
  selSIG=uSIG[i]#"ROS-010_1807-2307_POSP"
  mpr=ModelParams_Replacement %>% filter(SIG_REPLACEMENT==T,GENUS_CODE==MODs$GENUS_CODE[which(MODs$SIG==selSIG)[1]])
  SIGr=MODs %>% filter(SIG==selSIG)
  CurveRef=data.frame(SAx=y,SurvR=inv.logit(mpr$s.int+mpr$s.slp*y),GrowR=(mpr$pg.int+mpr$pg.slp*y),RecR=mpr$SecRecVal)
  RecRef=mpr$SecRecVal
  CurveSG=cbind(CurveRef,data.frame(SurvX=inv.logit(SIGr$s.int+SIGr$s.slp*y),GrowX=(SIGr$pg.int+SIGr$pg.slp*y)),RecX=SIGr$SecRecVal)
  
  ScurveP=CurveSG %>% ggplot() +
    geom_line(aes(x=SAx,y=SurvR),color="blue",size=2)+
    geom_line(aes(x=SAx,y=SurvX),color="gray",size=2)+
    theme_bw()
  GcurveP=CurveSG %>% ggplot() +
    geom_line(aes(x=SAx,y=GrowR),color="blue",size=2)+
    geom_line(aes(x=SAx,y=GrowX),color="gray",size=2)+
    geom_abline()+
    theme_bw()
  RcurveP=CurveSG %>% ggplot() +
    geom_line(aes(x=SAx,y=RecR),color="blue",size=2)+
    geom_line(aes(x=SAx,y=RecX),color="gray",size=2)+
    theme_bw()
  JoinP=ScurveP+GcurveP+RcurveP+
    plot_annotation(title=paste("SIG:",selSIG,"; L :",round(SIGr$Lambda,2),"; ",i,"of",length(uSIG),"; Ng:",SIGr$pg.N,"; Ns:",SIGr$s.N))
  print(JoinP)
 # readline(prompt="Press [enter] to continue")
}


pg.intM=matrix(rep(MODs$pg.int,length(y)),ncol=nrow(MODs),byrow=T)
pg.M=matrix(y,nrow=length(y),ncol=1)%*%matrix(MODs$pg.slp,nrow=1,ncol=nrow(MODs))+pg.intM
row.names(pg.M)=paste0("Y_pg_",1:length(y))

s.intM=matrix(rep(MODs$s.int,length(y)),ncol=nrow(MODs),byrow=T)
s.M=inv.logit(matrix(y,nrow=length(y),ncol=1)%*%matrix(MODs$s.slp,nrow=1,ncol=nrow(MODs))+s.intM)
row.names(s.M)=paste0("Y_s_",1:length(y))


CurvesDF=data.frame(MODs_pca[,c("REGION","SIG","GENUS_CODE","Lambda","REPLACE_B","SecRecVal","SecRecValSource")],t(pg.M),t(s.M)) 
CurvesDF$REPLACE_B=factor(CurvesDF$REPLACE_B,levels=c("ABOVE","BORDERLINE","BELOW"))
CurvesDF=CurvesDF %>%
  pivot_longer(cols=c(row.names(pg.M),row.names(s.M)),names_to = c("VitalRate","SAx_i"),names_prefix = "Y_",names_sep = "_",values_to = c("C")) %>% 
  mutate(SAx=y[as.numeric(SAx_i)]) %>% select(!SAx_i) %>% 
  pivot_wider(names_from = "VitalRate",values_from = "C")

CurvesDF %>% 
  ggplot(aes(x=SAx,y=s,color=Lambda,fill=Lambda,group=SIG))+
  geom_line()+
  scale_color_viridis_c()+
  scale_fill_viridis_c()+
  facet_grid(GENUS_CODE~REPLACE_B)


uSIG=unique(MODs$SIG)
uG=unique(MODs$GENUS_CODE)
#ROS-010_1807-2307_POSP
for (i in 1:length(uG)){
  selG=uG[i]#"ROS-010_1807-2307_POSP"
  mpr=ModelParams_Replacement %>% filter(SIG_REPLACEMENT==T,GENUS_CODE==selG)
  CurveRef=data.frame(SAx=y,SurvR=inv.logit(mpr$s.int+mpr$s.slp*y),GrowR=(mpr$pg.int+mpr$pg.slp*y),RecR=mpr$SecRecVal)

  ScurveP=CurvesDF %>% filter(GENUS_CODE==selG) %>% 
    ggplot() +
    geom_line(aes(x=SAx,y=s,color=REPLACE_B,group=SIG),size=1,alpha=0.5)+
    geom_line(aes(x=SAx,y=SurvR),color="blue",size=2,data=CurveRef)+
    scale_color_discrete(name="Pop. Replacement",direction=-1)+
    facet_grid(REPLACE_B~.)+
    theme_bw()
  GcurveP=CurvesDF %>% filter(GENUS_CODE==selG) %>% 
    ggplot() +
    geom_line(aes(x=SAx,y=pg,color=REPLACE_B,group=SIG),size=1,alpha=0.5)+
    geom_line(aes(x=SAx,y=GrowR),color="blue",size=2,data=CurveRef)+
    geom_abline()+
    scale_color_discrete(name="Pop. Replacement",direction=-1)+
    facet_grid(REPLACE_B~.)+
    theme_bw()
  RcurveP=CurvesDF %>% filter(GENUS_CODE==selG) %>% 
    ggplot() +
    geom_line(aes(x=SAx,y=SecRecVal,color=REPLACE_B,group=SIG),size=1,alpha=0.5)+
    geom_line(aes(x=SAx,y=RecR),color="blue",size=2,data=CurveRef)+
    scale_y_sqrt(limits=c(0,max(MODs$SecRecVal)))+
    scale_color_discrete(name="Pop. Replacement",direction=-1)+
    facet_grid(REPLACE_B~.)+
    theme_bw()
  JoinP=ScurveP+GcurveP+RcurveP+
    plot_layout(guides = "collect")+
    plot_annotation(title = paste0(selG," Vital Rates"))&theme(legend.position = "bottom")
  print(JoinP)
 # readline(prompt="Press [enter] to continue")
}

Elaslist
Name
save(list=c("MODs_pca","CurveRef","CurvesDF"),file=paste0("./Data/ModelData/",Name,"Model_LL_SIGCurves.rdata"))
     