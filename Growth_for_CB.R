library(tidyverse)
ct=load("./Data/ColonyTransitions/HA_Trans/Colony_Data_20210917_edited.rdata")
ct=ColonyLevel
ct = ct %>% mutate(LinExt=((EndingSize-StartingSize)/StartingPerim)/Interval_Years)  

ct %>% ggplot(aes(StartingSize,LinExt))+
  geom_point()+
  facet_wrap(Island~Genus_Code)+
  scale_x_log10()+
  #scale_y_log10()+ 
  theme_bw()

mnLE=ct %>%
  filter(TransitionTypeSimple=="GROWTH") %>%
  group_by(Island,Spec_Code) %>%
  reframe(Ncol=length(LinExt),
          mn_LinExt=mean(LinExt,na.rm=T),
          q50_LinExt=quantile(LinExt,na.rm=T,.5),
          q75_LinExt=quantile(LinExt,na.rm=T,.75),
          q90_LinExt=quantile(LinExt,na.rm=T,.9)) %>% 
  filter(Ncol>=6)

ct %>%   
  filter(TransitionTypeSimple%in%c("GROWTH","SHRINK")) %>%
  ggplot(aes(LinExt))+
  geom_histogram(aes(fill=TransitionTypeSimple))+
  facet_wrap(Island~Spec_Code,scale="free_y")+
  geom_vline(xintercept = 0,color="red")+
  geom_vline(aes(xintercept = q75_LinExt),color="purple",data=mnLE)+
  geom_vline(aes(xintercept = q90_LinExt),color="blue",data=mnLE)+
  scale_x_continuous(limits = c(-4,4))+
  #scale_x_log10()+
  #scale_y_log10()+ 
  theme_bw()


 View(mnLE)
 mnLE %>% ggplot(aes(Genus_Code,y=q90_LinExt))+geom_point()+facet_wrap("Island") 

 ct %>% ggplot(aes(StartingSize,EndingSize))+geom_point()+facet_wrap(Island~Genus_Code)+scale_x_log10()+scale_y_log10() 
 