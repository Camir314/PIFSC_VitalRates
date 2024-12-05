#Tom_Plotting
library(scales)
library(patchwork)
library(ggplot2)
library(tidyr)
library(dplyr)
library(lubridate)
library(fitdistrplus)
library(readr)
SmithsonVerkuilen2006=function(yi){
  n=length(yi)
  yo=(yi * (n-1) + 0.5) / n
  return(yo)
}
Inv.SmithsonVerkuilen2006=function(yo,n=length(yo)){
  yi=(n*yo - 0.5)/(n-1)
  return(yi)
}

#NEED THE A2D FUNCTION!!!!!


#############################
# Rec and Cover
#############################
#FIXED SITES

#Get to clean Area Surveyed
load("data/Colony_Data_20210622_edited.rdata")#Colony_Data_20210915_edited.rdata")
#HAWOCC010 and MAIOCC002 2016 date is wrong. Change it to correct date
#ColonyLevel$StartingDate[ColonyLevel$StartingDate == "2016-07-06"] <- "2016-08-11"
#ColonyLevel$StartingDate[ColonyLevel$StartingDate == "2016-12-19"] <- "2016-07-14"

Asurv=read.csv("./data/AreaSurveyed_N_Circrats.csv")
names(Asurv)[1]="Site"
Asurv$StartingDate=mdy(Asurv$Date)
Asurv_l=Asurv[,c("Site","StartingDate","POCS","POSP","MOSP")]%>%pivot_longer(cols=c("POSP","MOSP","POCS"),#,"dSE.POSP","dSE.MOSP","dSE.POCS"
                     names_to=c("Genus_Code"),values_to=c("Ncircrats"))
names(Asurv_l)[which(names(Asurv_l)=="StartingDate")]="EndingDate" #renaming StartingDate to EndingDate. Shouldn't line 34 just be EndingDate?

#Calculate Observed Recruitment per Surveyed Area
Rec.S_T=ColonyLevel %>% group_by(Site,StartingDate,EndingDate,Genus_Code,Interval_Years) %>% 
              summarize(Nrec=length(which(TransitionTypeSimple=="RECR"))) #count # of recruitment events for each site/year/genus in ColonyLevel
Rec.S_T.A=left_join(Rec.S_T,Asurv_l) #adds numb of circrats to the table created in last step
Rec.S_T.A$Ncircrats <- as.numeric(Rec.S_T.A$Ncircrats)
Rec.S_T.A$A.Surv.m2 =Rec.S_T.A$Ncircrats*0.5 #numb of circrats * 0.5m2 (size of circrat) to get area surveyed in m2
Rec.S_T.A$A.Surv.cm2 = Rec.S_T.A$Ncircrats*5000 #area surveyed cm2

#Get Area of adults by taxon in Surveyed Area
#Link to Nrecruits from the StartingSize summed area: Adults at beginning of interval generate 
Atax=ColonyLevel %>% group_by(Site,StartingDate,EndingDate,Genus_Code) %>% summarise(A.adult_start_cm2=sum(StartingSize)) #Area of adult cm2 = sum all corals in each genus for each year
Rec.S_T.A=left_join(Rec.S_T.A,Atax,by=c("Site","StartingDate","EndingDate","Genus_Code")) #matching by Site, EndDate, Genus

#Recruit N per adult area and annual rate
#Rec.S_T.A$Rec.m2=Rec.S_T.A$Nrec/Rec.S_T.A$A.Surv.m2 # numb recruits/area surveyed m2 = Rec.m2
Rec.S_T.A$Rec.Ad_cm2=Rec.S_T.A$Nrec/Rec.S_T.A$A.adult_start_cm2 # numb recruits/total adult area cm2
Rec.S_T.A$Rec.Ad_cm2_pYr=Rec.S_T.A$Rec.Ad_cm2/(Rec.S_T.A$Interval_Years) # numb recruits/total adult area cm2 annualized

#Best Estimates of Recruitment
SFM_Observed_Recruitment=Rec.S_T.A[,c("Site","StartingDate","EndingDate","Interval_Years","Genus_Code","Rec.Ad_cm2_pYr")]

#Drop "Hole Situation" Interval Years - i.e. find Intervals that span multiple time point with holes in the plot, and recognize that 
# recruit estimates coming from those colonies are likely not great, so let's drop them 
HoleSitch=ColonyLevel %>% group_by(Site,StartingDate,EndingDate,Interval_Years) %>% summarize(Ncol=length(StartingSize))
#If you have less than 30 colonies (of any TransitionType) for your interval year, you're a hole situation, Dropped.
HoleSitchUnder30=subset(HoleSitch,Ncol<30)
# plot(sort(HoleSitch$Ncol))
# sort(HoleSitch$Ncol)
#Given the list in HoleSitchUnder30, match and drop
drop_i=NULL
for(i in 1:nrow(HoleSitchUnder30)){
 thisdrop_i=which(SFM_Observed_Recruitment$Site==HoleSitchUnder30$Site[i]&
                  SFM_Observed_Recruitment$EndingDate==HoleSitchUnder30$EndingDate[i]&
                  SFM_Observed_Recruitment$Interval_Years==HoleSitchUnder30$Interval_Years[i])
 drop_i = unique(c(drop_i,thisdrop_i))
}
drop_i=c(drop_i, which(is.infinite(SFM_Observed_Recruitment$Rec.Ad_cm2_pYr)))
if(length(drop_i)>0){SFM_Observed_Recruitment=SFM_Observed_Recruitment[-drop_i,]}

dim(SFM_Observed_Recruitment)

#Now, Given Every unique site and interval in the ColonyLevel data, match (or not with best estimate from Rec)
SIG_withSFMRec=unique(ColonyLevel[,c("Site","StartingDate","EndingDate","Genus_Code")])
SIG_withSFMRec$GS_Rec_Adcm2_pYr=NA
for(i in 1:nrow(SIG_withSFMRec)){
  matchup=subset(SFM_Observed_Recruitment,
                 Site==SIG_withSFMRec$Site[i]&
                   EndingDate==SIG_withSFMRec$StartingDate[i]&
                   Genus_Code==SIG_withSFMRec$Genus_Code[i])
  if(nrow(matchup)==1){
    print("YAY")
    SIG_withSFMRec$GS_Rec_Adcm2_pYr[i]=matchup$Rec.Ad_cm2_pYr  
  }else if(nrow(matchup)>1){print("WHAT?!")}
}
table(is.na(SIG_withSFMRec$GS_Rec_Adcm2_pYr))
save(list = "SIG_withSFMRec",file = "./data/ObservedSFMRecruitment_by_Site_Interval_Genus.rdata")

# 
# ggplot(SFM_Observed_Recruitment,aes(x=Rec.Ad_cm2_pYr,fill=Genus_Code))+
#   geom_histogram()+#geom_violin(position="dodge")+
#   theme_bw()+facet_grid("Genus_Code")+
#   theme(axis.text.x = element_text(angle=90))+
#   scale_x_sqrt()
# 
# ggplot(SFM_Observed_Recruitment,aes(EndingDate,Rec.Ad_cm2_pYr))+
#   geom_vline(data=ColonyLevel,aes(xintercept=StartingDate),color="blue")+
#   geom_vline(data=ColonyLevel,aes(xintercept=EndingDate),color="red")+
#   geom_point(color="purple")+
#   facet_wrap(Site~Genus_Code)+
#   xlim(c(ymd(20120101),ymd(20200101)))
  

#Here's our final best estimate of larval recruitment per adult area per time.
#Let's (1) fits exponential distribution to the data, and (2) come up with a median of each taxon
MRec=ddply(SFM_Observed_Recruitment,.(Genus_Code),summarize,
           q05=quantile(Rec.Ad_cm2_pYr,.05),
           med=median(Rec.Ad_cm2_pYr),
           mean=mean(Rec.Ad_cm2_pYr),
           q95=quantile(Rec.Ad_cm2_pYr,.95))


ggplot(SFM_Observed_Recruitment,aes(Rec.Ad_cm2_pYr))+
  geom_histogram()+
  geom_vline(data=MRec,aes(xintercept = q05),color="grey",lty=2)+
  geom_vline(data=MRec,aes(xintercept = med),color="pink")+
  geom_vline(data=MRec,aes(xintercept = q95),color="grey",lty=2)+
  facet_grid(Genus_Code~.)


Rec.S_T.A=subset(Rec.S_T.A,!is.infinite(Rec.Ad_cm2_pYr))
par(mfrow=c(1,3))
descdist(subset(Rec.S_T.A,Genus_Code=="POCS")$Rec.Ad_cm2_pYr,boot = 1000)
descdist(subset(Rec.S_T.A,Genus_Code=="POSP")$Rec.Ad_cm2_pYr,boot = 1000)
descdist(subset(Rec.S_T.A,Genus_Code=="MOSP")$Rec.Ad_cm2_pYr,boot = 1000)

#All Looks roughly exponentially distributed
summary(Rec.S_T.A$Rec.Ad_cm2_pYr)
maxnorm=max(Rec.S_T.A$Rec.Ad_cm2_pYr)
Rec.S_T.A$Rec.Ad_cm2_pYr_b=SmithsonVerkuilen2006(Rec.S_T.A$Rec.Ad_cm2_pYr/maxnorm)
Npocs=nrow(subset(Rec.S_T.A,Genus_Code=="POCS"))
Nposp=nrow(subset(Rec.S_T.A,Genus_Code=="POSP"))
Nmosp=nrow(subset(Rec.S_T.A,Genus_Code=="MOSP"))
POCSb=fitdist(subset(Rec.S_T.A,Genus_Code=="POCS")$Rec.Ad_cm2_pYr_b,distr = "beta")
POSPb=fitdist(subset(Rec.S_T.A,Genus_Code=="POSP")$Rec.Ad_cm2_pYr_b,distr = "beta")
MOSPb=fitdist(subset(Rec.S_T.A,Genus_Code=="MOSP")$Rec.Ad_cm2_pYr_b,distr = "beta")
xs=seq(0,1,length.out = 1000)
y1=dbeta(x = xs,shape1 = POCSb$estimate[1],shape2 = POCSb$estimate[2])
q1=qbeta(p = c(.05,.2,.5,.8,.95),shape1 = POCSb$estimate[1],shape2 = POCSb$estimate[2])
m1t=maxnorm*Inv.SmithsonVerkuilen2006(1/(1+POCSb$estimate[2]/POCSb$estimate[1]),n=Npocs)
y2=dbeta(x = xs,shape1 = POSPb$estimate[1],shape2 = POSPb$estimate[2])
q2=qbeta(p = c(.05,.2,.5,.8,.95),shape1 = POSPb$estimate[1],shape2 = POSPb$estimate[2])
m2t=maxnorm*Inv.SmithsonVerkuilen2006(1/(1+POSPb$estimate[2]/POSPb$estimate[1]),n=Nposp)
y3=dbeta(x = xs,shape1 = MOSPb$estimate[1],shape2 = MOSPb$estimate[2])
q3=qbeta(p = c(.05,.2,.5,.8,.95),shape1 = MOSPb$estimate[1],shape2 = MOSPb$estimate[2])
m3t=maxnorm*Inv.SmithsonVerkuilen2006(1/(1+MOSPb$estimate[2]/MOSPb$estimate[1]),n=Nmosp)
xst=maxnorm*Inv.SmithsonVerkuilen2006(xs)
q1t=maxnorm*Inv.SmithsonVerkuilen2006(q1,n=Npocs)
q2t=maxnorm*Inv.SmithsonVerkuilen2006(q2,n=Nposp)
q3t=maxnorm*Inv.SmithsonVerkuilen2006(q3,n=Nmosp)

par(mfrow=c(3,1))
hist(subset(Rec.S_T.A,Genus_Code=="POCS")$Rec.Ad_cm2_pYr,30,main="POCS",xlab="Recruitment Rate (N/(cm2 adult * yr))",breaks=seq(0,.04,length.out=100))
points(xst,y1,type="l",col="grey25",lty=2)
abline(v=q1t[c(1,5)],col="grey75",lty=2)
abline(v=q1t[c(2,4)],col="grey25",lty=2)
abline(v=q1t[c(3)],col="black",lty=1)
abline(v=m1t,col="red",lty=1)
hist(subset(Rec.S_T.A,Genus_Code=="POSP")$Rec.Ad_cm2_pYr,30,main="POSP",xlab="Recruitment Rate (N/(cm2 adult * yr))",breaks=seq(0,.04,length.out=100))
points(xst,y2,type="l",col="grey25",lty=2)
abline(v=q2t[c(1,5)],col="grey75",lty=2)
abline(v=q2t[c(2,4)],col="grey25",lty=2)
abline(v=q2t[c(3)],col="black",lty=1)
abline(v=m2t,col="red",lty=1)
hist(subset(Rec.S_T.A,Genus_Code=="MOSP")$Rec.Ad_cm2_pYr,30,main="MOSP",xlab="Recruitment Rate (N/(cm2 adult * yr))",breaks=seq(0,.04,length.out=100))
points(xst,y3,type="l",col="grey25",lty=2)
abline(v=q3t[c(1,5)],col="grey75",lty=2)
abline(v=q3t[c(2,4)],col="grey25",lty=2)
abline(v=q3t[c(3)],col="black",lty=1)
abline(v=m3t,col="red",lty=1)

ms=c(m1t,m2t,m3t);names(ms)=c("POCS","POSP","MOSP")
ddply(Rec.S_T.A,.(Genus_Code),summarize,mn=mean(Rec.Ad_cm2_pYr))

#calculate mean/median juvenile density (number of recruits/total adult area)
Rec.S_T.A$Rec.Ad_cm2[is.infinite(Rec.S_T.A$Rec.Ad_cm2)] <- NA
mean(Rec.S_T.A$Rec.Ad_cm2, na.rm = TRUE) #overall mean
median(Rec.S_T.A$Rec.Ad_cm2, na.rm = TRUE) #overall median
FIXmeanMOSP <- mean(Rec.S_T.A$Rec.Ad_cm2[Rec.S_T.A$Genus_Code == "MOSP"]) 
FIXmeanPOSP <- mean(Rec.S_T.A$Rec.Ad_cm2[Rec.S_T.A$Genus_Code == "POSP"])
FIXmeanPOCS <- mean(Rec.S_T.A$Rec.Ad_cm2[Rec.S_T.A$Genus_Code == "POCS"], na.rm = TRUE)

hist(Rec.S_T.A$Rec.Ad_cm2)
FIX=ggplot(Rec.S_T.A,aes(Rec.Ad_cm2,fill=Genus_Code))+
  geom_histogram(binwidth = 0.005)+
  facet_grid("Genus_Code")+
  xlim(c(-.01,.15))+
  ggtitle("Fixed Sites")+
  xlab("Number recruits/total adult area")+
  theme(legend.position = "none")


#count # of recruitment events for each site/year/genus in ColonyLevel
recrui <- subset(ColonyLevel, TransitionTypeSimple == "RECR")
#need to edit this
numbevents=ddply(recrui,.(Site,StartingDate,Genus_Code,ln_ES),
              summarize,Nrec=length(which(TransitionTypeSimple=="RECR"))) 

#average recruit size = average all corals in each genus for each year
recrui$ln_ES[recrui$ln_ES == -Inf] <-NA
avgrecr_size =ddply(recrui,.(Site,StartingDate,Genus_Code),summarise,AvgEndSize=mean(ln_ES, na.rm = TRUE )) 

#plot recruit histogram
png(filename = "Figures/SizeFreqDist_Juveniles_Fixex.png")
ggplot(avgrecr_size, aes(x= AvgEndSize, fill = Genus_Code))+
  geom_histogram()+
  facet_wrap(~Genus_Code)+
  ggtitle("Recruit Size (Log Transformed)")+
  xlab("Log10 Size")+ ylab("Recruits")
dev.off()

################
#REA DATA
#Check Proportion of Recruits per Taxon - GGPLOT
RECs=subset(ColonyLevel,TransitionTypeSimple=="RECR")
bw=0.5
ggplot()+
  geom_histogram(data=subset(ColonyLevel,EndingSize>0),aes(A2D(EndingSize)),fill="grey",binwidth=bw)+
  xlim(c(0,10))+
  geom_histogram(data=RECs,aes(A2D(EndingSize)),fill="blue",binwidth=bw)+
  facet_grid(Genus_Code~.,scales = "free_y")+
  geom_vline(xintercept = 5)+
  theme_bw()

#Check Proportion of Recruits per Taxon - hist calculation
binwidths=bw#seq(.01,1,length.out=100)#.333333
uG=c("POCS","POSP","MOSP")
PropRec_g=cbind(expand.grid(Genus_Code=uG,BinWidth=binwidths),PropRec=NA)
#par(mfrow=c(3,1))
for(i_b in 1:length(binwidths)){
  binwidth=binwidths[i_b]
  N=5/binwidth
  for(i_G in 1:length(uG)){
    hcl=hist(A2D(subset(ColonyLevel,EndingSize>0&Genus_Code==uG[i_G])$EndingSize),breaks=seq(0,1000,by=binwidth),plot = F)
    hrc=hist(A2D(subset(RECs,Genus_Code==uG[i_G])$EndingSize),breaks=seq(0,1000,by=binwidth),plot = F)
    PropRec_v=(hrc$counts[1:(N+1)]/hcl$counts[1:(N+1)])
    out_i=which(PropRec_g$Genus_Code==uG[i_G]&PropRec_g$BinWidth==binwidth)
    PropRec_g$PropRec[out_i]=mean(PropRec_v[1:(N+1)],na.rm=T)
    # plot(hcl$breaks[1:(N+1)],PropRec_v,type="b",ylim=c(0,1),main=uG[i_G])
    # abline(h=mean(propREC,na.rm=T))
  }
}

PropRecMn=PropRec_g %>% 
  group_by(Genus_Code) %>% 
  summarise(mnPR=mean(PropRec))

#Final 0-5 cm diam. Proportion Recruits in the Juvenile size classes
PropRecMn


#Juv Data, Cover Data
juv <- read.csv("./data/BenthicREA_sectordata_GENUS.csv")
juv_site=read.csv("./data/JuvProject_temporal_SITE.csv")
cov=read.csv("./data/BenthicCover_2010-2019_Tier3_SECTOR.csv")
cov_site=read.csv("./data/BenthicCover_2010-2020_Tier3_SITE.csv")
jmhi=subset(juv,REGION=="MHI"&GENUS_CODE%in%c("POCS","MOSP","POSP"))
jmhi_site=subset(juv_site,REGION=="MHI"&GENUS_CODE%in%c("POCS","MOSP","POSP"))
cmhi=subset(cov,Mean.REGION=="MHI")
cmhi_site=subset(cov_site,REGION=="MHI")

#Get Sector Level Cover, POSP, MOSP, POCS # 5 cm diameter.
TaxSec_i=which(substr(names(cmhi),1,5)=="Mean."&!substr(names(cmhi),6,50)%in%c("REGION","ISLAND","ANALYSIS_SEC","ANALYSIS_YEAR","TOT_AREA_WT","N"))
names(cmhi)[TaxSec_i]
rowSums(cmhi[,TaxSec_i])
cmhi$Mean.MOSP=rowSums(cmhi[,c("Mean.MOEN","Mean.MOBR","Mean.MOFO")])#"Mean.MONE",
cmhi$PooledSE.MOSP=sqrt(rowSums((sqrt(cmhi$Mean.N)*cmhi[,c("PooledSE.MOEN","PooledSE.MOBR","PooledSE.MOFO")])^2))/sqrt(cmhi$Mean.N)#"Mean.MONE",
cmhi$Mean.POSP=rowSums(cmhi[,c("Mean.POMA","Mean.POEN")])
cmhi$PooledSE.POSP=sqrt(rowSums((sqrt(cmhi$Mean.N)*cmhi[,c("PooledSE.POMA","PooledSE.POEN")])^2))/sqrt(cmhi$Mean.N)
metacol_sec=c("Mean.REGION","Mean.ISLAND","Mean.ANALYSIS_SEC","Mean.ANALYSIS_YEAR","Mean.N")
cov_sec_w=cmhi[,c(metacol_sec,"Mean.POSP","Mean.MOSP","Mean.POCS","PooledSE.POSP","PooledSE.MOSP","PooledSE.POCS")]
names(cov_sec_w)=substr(names(cov_sec_w),6,99)
cov_long=cov_sec_w[,names(cov_sec_w)[1:8]]%>% 
  pivot_longer(cols=c("POSP","MOSP","POCS"),#,"dSE.POSP","dSE.MOSP","dSE.POCS"
               names_to=c("GENUS_CODE"),values_to=c("cover.mn"))
# cov_sd=cov_sec_w[,c(names(cov_sec_w)[1:5],c("dSE.POSP","dSE.MOSP","dSE.POCS"))] %>% 
#   pivot_longer(cols=c("dSE.POSP","dSE.MOSP","dSE.POCS"),
#                names_to=c("dSE.Genus"),values_to=c("cover.se"))
# cov_sec=left_join(cov_long,cov_sd)
cov_sec=cov_long#select(cov_sec,-dSE.Genus)

names(jmhi)[which(names(jmhi)=="Sector")]="ANALYSIS_SEC"
juv_sec=jmhi[,c("REGION","ISLAND","ANALYSIS_SEC","ANALYSIS_YEAR","GENUS_CODE","n","Mean_JuvColDen","SE_JuvColDen")]
names(juv_sec)[names(juv_sec)=="n"]="N.juv"
names(cov_sec)[names(cov_sec)=="N"]="N.cov"
cov_sec$ANALYSIS_YEAR[cov_sec$ANALYSIS_YEAR=="2013-15"]="2013"
cov_sec=subset(cov_sec,ANALYSIS_YEAR%in%c("2013","2016","2019"))
cov_sec$ANALYSIS_YEAR=as.numeric(as.vector(cov_sec$ANALYSIS_YEAR))
sec=na.omit(left_join(juv_sec,cov_sec,by=c("REGION","ISLAND","ANALYSIS_SEC","ANALYSIS_YEAR","GENUS_CODE")))

sec$coverP.mn=sec$cover.mn/100 #as prop
#sec$coverP.se=sec$cover.se/100
#sec$coverP.var=(sqrt(sec$N.cov)*sec$coverP.se)^2
sec$juv.var=(sqrt(sec$N.juv)*sec$SE_JuvColDen)^2
sec$JuvColDen.mn_cm2=sec$Mean_JuvColDen/10000 #convert to cm2
sec$juv_mod.mn.pcm2=sec$JuvColDen.mn_cm2/(sec$coverP.mn) #mean juv den/mean cover
#sec$juv_mod.var.pcm2=sec$juv_mod.mn.pcm2*sqrt((sec$juv.var/sec$Mean_JuvColDen)^2+(sec$coverP.var/sec$coverP.mn)^2)
#sec$juv_mod.se.pcm2=sqrt(sec$juv_mod.var)/sqrt(sec$N.juv) #Juvenile density (number of recruits/adult cm2 area)


REA=ggplot(sec,aes(x = juv_mod.mn.pcm2,fill=GENUS_CODE))+
  geom_histogram(binwidth = .005)+
  #geom_errorbar()+
  facet_grid(GENUS_CODE~.)+
  xlim(c(-.01,.15))+
  xlab("Juvenile Density per Area Coral Cover")+
  ggtitle("REA data")+
  theme(axis.title.y = element_blank())

smallguys=subset(ColonyLevel,log10_ES<(log10(2.5^2*pi))&TransitionType%in%c("MORT","RECR","GROWTH","SHRINK","FISSION","FUSION"))
t=table(smallguys$Genus_Code,smallguys$TransitionType)
PercentREC=t[,"RECR"]/rowSums(t)

sec$juv_mod.mn.pcm2_prorate=sec$juv_mod.mn.pcm2*PercentREC[sec$GENUS_CODE]
  
REApro=ggplot(sec,aes(x = juv_mod.mn.pcm2_prorate,fill=GENUS_CODE))+
  geom_histogram(binwidth = .005)+
  #geom_errorbar()+
  facet_grid(GENUS_CODE~.)+
  xlim(c(-.01,.15))+
  xlab("Juvenile Density per Adult cm^2")+
  ggtitle("REA data - PRORATED By Prop. Recruit")

FIX+REA+REApro

sec %>% group_by(GENUS_CODE) %>% summarize(med=median(juv_mod.mn.pcm2))

png(filename = "Figures/ObservedJuvDen.png")
FIX+REA
dev.off()


#overall mean
mean(sec$juv_mod.mn.pcm2, na.rm = TRUE)
median(sec$juv_mod.mn.pcm2, na.rm = TRUE)
#genus
meanMOSP <- mean(sec$juv_mod.mn.pcm2[sec$GENUS_CODE == "MOSP"])
meanPOSP <- mean(sec$juv_mod.mn.pcm2[sec$GENUS_CODE == "POSP"])
meanPOCS <- mean(sec$juv_mod.mn.pcm2[sec$GENUS_CODE == "POCS"])

#recruit size histogram







ggplot(sec,aes(x = ANALYSIS_SEC,y =juv_mod.mn.pcm2,ymin =juv_mod.mn.pcm2-juv_mod.se.pcm2,ymax =juv_mod.mn.pcm2+juv_mod.se.pcm2,fill=ISLAND))+
  geom_point()+
  #geom_errorbar()+
  facet_grid(.~GENUS_CODE)+
  theme(axis.text.x = element_text(angle=90))+ylab("Juvenile Density per Adult cm^2 area")+ylim(c(0,.3))

ggplot(sec,aes(mean = juv_mod.mn.pcm2, sd = sqrt(juv_mod.var.pcm2),fill=ISLAND))+
  stat_function(fun = dnorm, geom = "area", 
                          fill = "orange", alpha = 0.25)+
  facet_grid("GENUS_CODE")

#scaled juvcolden
sec$Genus=GeneraLU[sec$GENUS_CODE]
sec$Genus=factor(sec$Genus,levels=c("Porites sp.","Montipora sp.","Pocillopora sp." ))
ggplot(sec,aes(juv_mod))+geom_histogram()+facet_grid("GENUS_CODE")+scale_x_log10()

#raw juvcolden, sector
ggplot(jmhi,aes(x=GENUS_CODE,y=Mean_JuvColDen))+geom_violin()+
  #facet_grid("GENUS_CODE")+
  scale_y_log10()
#raw juvcolden, site
ggplot(jmhi_site,aes(JuvColDen ))+geom_histogram()+facet_grid(GENUS_CODE~ISLAND)+scale_x_log10()


#Get Site Level Cover, POSP, MOSP, POCS
TaxSite_i=22:126#which(substr(names(cmhi_site),1,5)=="Mean."&substr(names(cmhi),6,16)!="TOT_AREA_WT")
names(cmhi_site)[TaxSite_i]
rowSums(cmhi_site[,TaxSite_i])
cmhi_site$MOSP=rowSums(cmhi_site[,c("MOEN","MOBR","MOFO")])#"Mean.MONE",
cmhi_site$POSP=rowSums(cmhi_site[,c("POMA","POEN")])
metacol_site=c("SITEVISITID","REGION","ISLAND","SEC_NAME","SITE","LATITUDE","LONGITUDE","ANALYSIS_YEAR","OBS_YEAR","DATE_","DEPTH_BIN")
site_cov=cmhi_site[,c(metacol_site,"POSP","MOSP","POCS")]

apply(cmhi[,5:ncol(cmhi)],2,mean)

head(cmhi)

ddply(jmhi,.(GENUS_CODE),summarize,
      q05=quantile(113*JuvColDen,.05),
      q25=quantile(113*JuvColDen,.25),
      q50=quantile(113*JuvColDen,.50),
      q75=quantile(113*JuvColDen,.75),
      q95=quantile(113*JuvColDen,.95))

ggplot(jmhi,aes(x=JuvColDen))+
  geom_histogram()+
  scale_x_log10()+
  facet_grid(GENUS_CODE~OBS_YEAR)+
  geom_vline(xintercept = 35.75*.1,color="red")+
  geom_vline(xintercept = 30.75*.1,color="red")


#############################
# Data Transitions
#############################
ColonyLevel$Genus=factor(ColonyLevel$Genus,levels=c("Porites sp.","Montipora sp.","Pocillopora sp." ))
ColonyLevel$TransitionTypeSimple=factor(ColonyLevel$TransitionTypeSimple,levels=c("GROWTH","SHRINK","RECR","MORT"))

zeq=0.01
TaxaTrans=ggplot(ColonyLevel,aes(x=StartingSize,y=AnnualEndingSize_E,
                       color=TransitionTypeSimple,
                       fill=TransitionTypeSimple,shape=TransitionTypeSimple))+
  geom_abline(intercept=0,slope=1)+
  geom_vline(xintercept = 1,lty=3,color="gray50")+
  geom_hline(yintercept = 1,lty=3,color="gray50")+
  geom_vline(xintercept = zeq,lty=1,color="black")+
  geom_hline(yintercept = zeq,lty=1,color="black")+
  geom_jitter(aes(x=StartingSize+zeq/2,y=EndingSize+zeq/2),height=0,width=10*zeq,
              data=subset(ColonyLevel,TransitionTypeSimple%in%c("RECR")),size=.5)+
  geom_point(size=1,data=subset(ColonyLevel,TransitionTypeSimple%in%c("GROWTH","SHRINK")))+
  geom_jitter(aes(x=StartingSize+zeq/2,y=EndingSize+zeq/2),height=10*zeq,width=0,
              data=subset(ColonyLevel,TransitionTypeSimple%in%c("MORT")),size=.5)+
  #geom_point(y=.001,aes(x=StartingSize),data=subset(ColonyLevel,TransitionTypeSimple=="MORT"),size=.5)+
  facet_grid(.~Genus)+
  scale_x_log10(name="Colony Area\nat T (cm^2)",labels=label_comma(drop0trailing =T),
                breaks=c(.1,1,10,100,1000,10000),limits=c(zeq/3,50000))+
  scale_y_log10(name="Colony Area\n at T + 1 Year (cm^2)",labels=label_comma(drop0trailing =T),
                breaks=c(.1,1,10,100,1000,10000),limits=c(zeq/3,50000))+
  scale_shape_manual(name="Transition Type",values=c(24,25,3,4),
                     breaks=c("GROWTH","SHRINK","RECR","MORT"),
                     labels=c("Growth","Partial Mortality","Recruitment","Mortality"))+
  scale_color_manual(name="Transition Type",values=c("darkblue","darkcyan","darkgreen","darkred"),
                     breaks=c("GROWTH","SHRINK","RECR","MORT"),
                     labels=c("Growth","Partial Mortality","Recruitment","Mortality"))+
  scale_fill_manual(name="Transition Type",values=c("darkblue","darkcyan","darkgreen","darkred"),
                    breaks=c("GROWTH","SHRINK","RECR","MORT"),
                    labels=c("Growth","Partial Mortality","Recruitment","Mortality"))+
  coord_equal()+theme_bw()+ggtitle("Colony-Level Transitions: All Sites, Intervals")+
  theme(legend.position = "bottom",axis.text=element_text(size=8))
TaxaTrans
sc=.75
ggsave(filename = "./figs/TaxaTransitions_All.png",TaxaTrans,width=sc*16,height=sc*9)


CL_POSP=subset(ColonyLevel,Genus_Code=="POSP")
SiteTrans=ggplot(CL_POSP,aes(x=StartingSize,y=AnnualEndingSize_E,
                                 color=TransitionTypeSimple,
                                 fill=TransitionTypeSimple,shape=TransitionTypeSimple))+
  geom_abline(intercept=0,slope=1)+
  geom_vline(xintercept = 1,lty=3,color="gray50")+
  geom_hline(yintercept = 1,lty=3,color="gray50")+
  geom_vline(xintercept = zeq,lty=1,color="black")+
  geom_hline(yintercept = zeq,lty=1,color="black")+
  geom_point(size=1,data=subset(CL_POSP,TransitionTypeSimple%in%c("GROWTH","SHRINK")))+
  geom_jitter(aes(x=StartingSize+zeq/2,y=EndingSize+zeq/2),height=10*zeq,width=0,
              data=subset(CL_POSP,TransitionTypeSimple%in%c("MORT")),size=.5)+
  geom_jitter(aes(x=StartingSize+zeq/2,y=EndingSize+zeq/2),height=0,width=10*zeq,
              data=subset(CL_POSP,TransitionTypeSimple%in%c("RECR")),size=.5)+
  facet_wrap(facets = c("Site"),nrow = 3,ncol = 5)+
  scale_x_log10(name="Colony Area\nat T (cm^2)",labels=label_comma(drop0trailing =T),
                breaks=c(.1,1,10,100,1000,10000),limits=c(zeq/3,50000))+
  scale_y_log10(name="Colony Area\n at T + 1 Year (cm^2)",labels=label_comma(drop0trailing =T),
                breaks=c(.1,1,10,100,1000,10000),limits=c(zeq/3,50000))+
  scale_shape_manual(name="Transition Type",values=c(24,25,3,4),
                     breaks=c("GROWTH","SHRINK","RECR","MORT"),
                     labels=c("Growth","Partial Mortality","Recruitment","Mortality"))+
  scale_color_manual(name="Transition Type",values=c("darkblue","darkcyan","darkgreen","darkred"),
                     breaks=c("GROWTH","SHRINK","RECR","MORT"),
                     labels=c("Growth","Partial Mortality","Recruitment","Mortality"))+
  scale_fill_manual(name="Transition Type",values=c("darkblue","darkcyan","darkgreen","darkred"),
                    breaks=c("GROWTH","SHRINK","RECR","MORT"),
                    labels=c("Growth","Partial Mortality","Recruitment","Mortality"))+
  coord_equal()+theme_bw()+ggtitle("Colony-Level Transitions: Porites, Each Site, All Intervals")+
  theme(legend.position = "bottom",axis.text=element_text(size=8))
sc=.75
ggsave(filename = "./figs/Porites_SiteTransitions_All.png",SiteTrans,width=sc*16,height=sc*9)


CL_POSP_SITE=subset(ColonyLevel,Genus_Code=="POSP"&Site=="MAI_SIO_K01")
SITrans=ggplot(CL_POSP_SITE,aes(x=StartingSize,y=AnnualEndingSize_E,
                             color=TransitionTypeSimple,
                             fill=TransitionTypeSimple,shape=TransitionTypeSimple))+
  geom_abline(intercept=0,slope=1)+
  geom_vline(xintercept = 1,lty=3,color="gray50")+
  geom_hline(yintercept = 1,lty=3,color="gray50")+
  geom_vline(xintercept = zeq,lty=1,color="black")+
  geom_hline(yintercept = zeq,lty=1,color="black")+
  geom_point(size=1,data=subset(CL_POSP_SITE,TransitionTypeSimple%in%c("GROWTH","SHRINK")))+
  geom_jitter(aes(x=StartingSize+zeq/2,y=EndingSize+zeq/2),height=10*zeq,width=0,
              data=subset(CL_POSP_SITE,TransitionTypeSimple%in%c("MORT")),size=.5)+
  geom_jitter(aes(x=StartingSize+zeq/2,y=EndingSize+zeq/2),height=0,width=10*zeq,
              data=subset(CL_POSP_SITE,TransitionTypeSimple%in%c("RECR")),size=.5)+
  facet_wrap(facets = c("Interval"))+
  scale_x_log10(name="Colony Area\nat T (cm^2)",labels=label_comma(drop0trailing =T),
                breaks=c(.1,1,10,100,1000,10000),limits=c(zeq/3,50000))+
  scale_y_log10(name="Colony Area\n at T + 1 Year (cm^2)",labels=label_comma(drop0trailing =T),
                breaks=c(.1,1,10,100,1000,10000),limits=c(zeq/3,50000))+
  scale_shape_manual(name="Transition Type",values=c(24,25,3,4),
                     breaks=c("GROWTH","SHRINK","RECR","MORT"),
                     labels=c("Growth","Partial Mortality","Recruitment","Mortality"))+
  scale_color_manual(name="Transition Type",values=c("darkblue","darkcyan","darkgreen","darkred"),
                     breaks=c("GROWTH","SHRINK","RECR","MORT"),
                     labels=c("Growth","Partial Mortality","Recruitment","Mortality"))+
  scale_fill_manual(name="Transition Type",values=c("darkblue","darkcyan","darkgreen","darkred"),
                    breaks=c("GROWTH","SHRINK","RECR","MORT"),
                    labels=c("Growth","Partial Mortality","Recruitment","Mortality"))+
  coord_equal()+theme_bw()+ggtitle("Colony-Level Transitions: Porites, Maui, Olowalu #3, All Intervals")+
  theme(legend.position = "bottom",axis.text=element_text(size=8))
sc=.7
ggsave(filename = "./figs/Porites_OL3.png",SITrans,width=sc*16,height=sc*9)

#RecVal All Taxa
recvalpc <- 0.32599#35.53
recvalpp <- 2.229822#30.61
recvalmp <- 0.1669464#31.74 

recdf=data.frame(Taxon=c("Porites sp.","Montipora sp.","Pocillopora sp." ),rec=c(recvalpp,recvalmp,recvalpc))
recdf$Taxon=factor(recdf$Taxon,levels=c("Porites sp.","Montipora sp.","Pocillopora sp." ))
recplot=ggplot(recdf,aes(x=Taxon,y=rec))+
  geom_col(width=.5,fill="skyblue")+
  ylab("N. Recruits\n (relative to adult area)")+
  ggtitle("Recruitment Tuned to Population Replacement (Lambda = 1)")+theme_classic()
sc=.7
ggsave(filename = "./figs/RecruitComparison_tunedLam1.png",recplot,width=sc*16,height=sc*9)



