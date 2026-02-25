library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

CBdb=read.csv("C:/Users/Thomas.Oliver/Downloads/prod_dbase_ha.csv")
head(CBdb)

DENdb=read.csv("C:/Users/Thomas.Oliver/Downloads/export_20260209/data_20260209.csv")
head(DENdb)

Nsp=sort(table(DENdb$specie_name))
Nsp10=names(Nsp[Nsp>=10])
DENdb=DENdb %>% filter(trait_name=="Skeletal density") %>%
  mutate(sden=as.numeric(value),genus=substr(specie_name,1,regexpr(pattern = " ",text = specie_name)))
genmed=DENdb %>% group_by(genus) %>% reframe(sden_med=median(sden)) %>% arrange(desc(sden_med))
DENdb=DENdb %>% 
  mutate(genus=factor(genus,levels = genmed$genus))

DENdb.=DENdb %>% filter(specie_name%in%Nsp10,!is.na(latitude))


DENdb %>% ggplot(aes(x=specie_name,y=sden))+geom_boxplot()+theme(axis.text.x = element_text(angle=90))+facet_wrap("genus",scales="free_x",drop=T)
DENdb %>% ggplot(aes(x=genus,y=sden))+geom_boxplot()+theme(axis.text.x = element_text(angle=90))
DENdb. %>% ggplot(aes(x=sden))+geom_histogram()+facet_wrap("specie_name")
DENdb. %>% ggplot(aes(x=specie_name,y=sden))+geom_point()

world <- ne_countries(scale = "medium", returnclass = "sf")
DENdb.sp=st_as_sf(DENdb.,coords = c("longitude","latitude"),crs=st_crs(world))
DENdbsp=st_as_sf(DENdb.,coords = c("longitude","latitude"),crs=st_crs(world))


DENdbsp %>% ggplot()+
  geom_sf(fill = "light gray",data = world) +
  geom_sf(aes(color=specie_name,size=sden))+
  scale_size_area()+
  theme(panel.background = element_rect(fill = "azure"))



sort(table(DENdb$specie_name))


mod1=lm(sden~specie_name+latitude,data=DENdb.)

summary(mod1)


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
 