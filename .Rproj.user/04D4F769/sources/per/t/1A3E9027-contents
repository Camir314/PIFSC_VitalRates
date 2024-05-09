#### MARAMP 2022 Vital Rates QC, initial dataframe creations and descriptive stats
### Corinne Amir 
### July 2023

library(dplyr)
library(anytime)
library(ggplot2)
library(ggmap)
library(viridis)
library(ggspatial)
library(ggrepel)
library(stringr)
library(sp)


setwd('C:/Users/Corinne.Amir/Documents/GitHub/PIFSC_VitalRates/CSV files') # Github repo

# raw <- read.csv("MARAMP22_VitalRates_07-24-2023.csv")
raw <- read.csv("ASRAMP23_VitalRates_05-07-2024.csv")
ll <- read.csv("VitalRates_LatLong.csv")
effort <- read.csv("VitalRates_SurveyEffort.csv")



#### QC Data ####
## (OPTIONAL) Remove superfluous columns:
colnames(raw)
raw <- raw %>% select(-c(OID_, TL_SurfA,QC_Check))


## Look for potential issues in the data:
lapply(raw, unique)



# ## Remove or alter rows based on TL_Note:
# raw <- raw %>% filter(TL_Note != "Out of bounds in 2022") # remove colonies that could not be tracked into 2022
# 
# 
# ## QC Morph code:
# a <- raw[raw$Morph_Code != "MD" & raw$Morph_Code != "EM" & raw$Morph_Code != "BR" & raw$Morph_Code != "TB" 
#         & raw$Morph_Code != "PL" & raw$Morph_Code != "KN",]
# 
# raw <- raw %>% mutate(Morph_Code = case_when(Morph_Code == "RM" ~ "EM",
#                                              Morph_Code == "EM" ~ "EM",
#                                              Morph_Code == "TB" ~ "TB",
#                                              Morph_Code == "BR" ~ "BR",
#                                              Morph_Code == "KN" ~ "KN",
#                                              Morph_Code == "PL" ~ "PL",
#                                              Morph_Code == "MD" ~ "MD"))


## Check if TimePt is labelled correctly:
a <- raw %>% group_by(Site, TL_Date, TimePt) %>% summarise(sum(TimePt))
  # HOW-005-2017/2018 = Jan 1
  # TUT-019-2018 = Jan 1
# raw$TL_Date <- as.factor(raw$TL_Date)   
# raw <- raw %>% mutate(TL_Date = recode(TL_Date,"1/1/2017" = "4/11/2017"))
# raw$TL_Date <- raw$TL_Date[raw$Site == "OCC-HOW-005" & raw$TL_Date == "1/1/2018"] = "6/8/2018"
# raw$TL_Date <- raw$TL_Date[which(raw$Site == "OCC-TUT-017" & raw$TL_Date == "1/1/2018")] <- "7/6/2018"

raw$TL_Date <- anydate(raw$TL_Date) # Change date format


# Add leading zeros to genet and colony code

raw$TL_id <- str_pad(raw$TL_id, 3, pad = "0")
raw$TL_Genet <- str_pad(raw$TL_Genet, 3, pad = "0")

#### Create additional columns #####

vr <- raw

# Turn TimePt into Year
vr$Year <- str_sub(vr$TL_Date,1,4)
vr %>% group_by(Site, TL_Date, Year) %>% sumvrse() # QC check



## Create unique name for all genets: 
vr$Genet_full <- paste(vr$Site,  vr$TL_Genet, vr$Year, sep = "_")
vr[11,] # double check

## Create unique name for all patches: 
vr$Patch_full <- paste(vr$Site, vr$TL_id, vr$Year,  sep = "_")
vr[11,] # double check


## Roll TL_Class up to genus (consider separating PGRA)
vr <- vr %>% mutate(Genus = case_when(TL_Class == "PMEA" ~ "POCS",
                                          TL_Class == "PVER" ~ "POCS",
                                          TL_Class == "AGLO" ~ "ACSP",
                                          TL_Class == "PGRA" ~ "POCS",
                                          TL_Class == "PLOB" ~ "POSP",
                                          TL_Class == "PLUT" ~ "POSP",
                                          TL_Class == "POCS" ~ "POCS",
                                          TL_Class == "MOSP" ~ "MOSP",
                                          TL_Class == "ACSP" ~ "ACSP",
                                          TL_Class == "POSP" ~ "POSP"))
lapply(vr, unique) # double check


## Add Lat and Long
ll <- ll %>% filter(Region != "vr") %>% rename(Site = ESD.Site.Name) %>% select(-Year)

vr <- left_join(vr, ll)


## Add m2 surveyed (collected from tracking spreadsheet)

vr <- left_join(vr, effort)


## Add area:perimeter ratio

vr$area_perim <- vr$Shape_Area/vr$Shape_Leng

head(vr)

## Add region code

# vr <- vr %>% mutate(Region = case_when(Island == "MAU" ~ "North",
#                                            Island == "PAG" ~ "North",
#                                            Island == "ASC" ~ "North",
#                                            Island == "SAI" ~ "South",
#                                            Island == "GUA" ~ "South"))


#### Descriptive Tables using patch data ####

a <- vr %>% group_by(Site, Annotator, TimePt) %>% summarise(n = n()) # Site-Years with the most and least patches
# OCC-SAI-009_2022: 352
# OCC-MAU-002_2017: 385
# OCC-GUA-015_2022: 25
# OCC-PAG-013_2017: 59


vr %>% group_by(Annotator) %>% summarise(n = n()) # Total patches surveyed by each annotator 
# JC: 1475 / 7 = 211 patches per Site-Year
# CA: 1127 / 9 = 125
# MSL: 539 / 4 = 135
# IGB: 428 / 4 = 107


vr %>% group_by(Genus) %>% summarise(n = n()) # How many total patches per genus
# POSP: 2702
# ACSP: 402
# POCS: 383
# MOSP: 82

vr %>% group_by(Genus, TimePt) %>% summarise(n = n()) # Patches per species over time (account survey effort)
# POSP: 1083 -> 1356 -> 263
# ACSP: 233 -> 145 -> 24
# POCS: 183 -> 90 -> 110
# MOSP: 75 -> 6 -> 1


a <- vr %>% group_by(Site, Genus, TimePt) %>% summarise(n = log(mean(TL_Area))) # Mean patch size

a <- vr %>% group_by(Genus, Year, Island) %>% summarise(n = n())

a <- vr %>% group_by(Year,Genus) %>% summarise(n = n())

#### Consolidate into colony dataframe ####

# Create colony dataframe:
vr_col <- vr %>% dplyr::select(Genet_full, TL_Area, TL_Perim, Shape_Leng, Shape_Area, area_perim) 
vr_meta <- vr %>% dplyr::select(Genet_full, Island, Region, Site, TimePt, Year, TL_Date, Latitude,Longitude,
                             Genus, TL_Class, TL_Genet, Quadrat, Effort) %>%
                      distinct()
vr_col <- aggregate(.~Genet_full, data = vr_col, sum)
vr_col <- left_join(vr_meta,vr_col)

              #group_by(Genet_full) %>% sumvrse(nrow(patches)) # Add in patch count 

# Total genets per Site-Year (ONLY colonies that survived 2+ time points):
vr_col$Site_Genet <- paste(vr_col$Site,vr_col$TL_Genet,  sep = "_") # use this column to filter


# Add column for number of patches
a <- vr %>% group_by(Site, Year,TL_Genet) %>% summarise(nPatches = n())

vr_col <- left_join(vr_col, a)

#### Descriptive Tables using colony data ####
# Total genets per Site-Year:
a <- vr_col %>% group_by(Site, TimePt, Genus) %>% 
        summarise(nColonies = n())

a <- vr_col %>% filter(TimePt !=0) %>% # look for recruits
                  group_by(Site_Genet) %>% 
                  filter(n()==1) %>% 
                  group_by(Site, Year,Genus) %>%
                  summarise(nColonies = n()) 

a <- vr_col %>% filter(TimePt !=2) %>% # look for growth/shrinkage/fission/fusion events
                  group_by(Site_Genet) %>% 
                  filter(n()>1) %>% 
                  group_by(Site, Genus) %>%
                  summarise(nColonies = n()) %>%
                  mutate(actual = nColonies/2)

a <- vr_col %>% filter(Island == "MAU" & TimePt !=0) %>% # find #transitions for additional Maug time point 
                  droplevels() %>%
                  group_by(Site_Genet) %>% 
                  filter(n()>1) %>% 
                  group_by(Site, Genus) %>%
                  summarise(nColonies = n()) %>%
                  mutate(actual = nColonies/2)

# Total area annotated per Site-Year:
b <- vr_col %>% group_by(Site, TimePt) %>% 
        summarise(tot_area = sum(Shape_Area))




# Deal with dispute in number of colonies in mult time points:
# a <- vr_change %>% group_by(Site, Genus)%>%
#   summarise(nColonies = n())
# a <- vr_col %>% filter(TimePt !=2) %>% # adding TimePt filter shows just 2014-2017 transition for Maug Sites
#   group_by(Site_Genet) %>% 
#   filter(n()>1) %>% 
#   group_by(Site, Genus) %>%
#   summarise(nColonies = n()) %>%
#   mutate(actual = nColonies/2)
 
#### Preliminary Plots ####

# Size frequency distribution (patches)
ggplot(data = mari %>% filter(Site =="OCC-PAG-006"), 
       aes(x=log(TL_Area))) + # Total
  geom_histogram() +
  facet_wrap(vars(Genus,Year), scales = "free_y", nrow = 3) 


# Kernel Density Plots (patches)
# mline <- mari %>% group_by(Site, Genus, Year) %>% summarise(Mean = log(mean(TL_Area))) # mean value by site
mline <- mari %>% group_by(Region, Genus, Year) %>% summarise(Mean = log(mean(Shape_Area))) # mean value by region


aa <-ggplot(data = mari %>% filter(Region == "South" & Year != "2014") %>% droplevels(), # By Site (plug and chug)
       aes(x=log(Shape_Area), group = Year, fill = Year)) +
  # geom_histogram(aes(y = ..density..),alpha = 0.4, binwidth = .4) + overlay histogram (it isnt raw counts, its proportion)
  geom_density() +
  # scale_fill_manual(values = alpha(c("#F1A340", "#998EC3"),0.75)) +
  scale_fill_manual(values = alpha(c("#F1A340", "#998EC3","aquamarine1"),0.75)) + # Maug sites
  facet_wrap(vars(Genus), nrow = 1, drop = F) + 
  geom_vline(data = mline %>% filter(Region == "South"& Year != "2014") ,
             aes(xintercept = Mean, color = Year),
             size = 1.1 ) +
  # scale_color_discrete(c("chocolate", "darkcyan", "chartreuse3")) +
  ggtitle("South/Populated") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())



# Other stuff I wanna plot
mari %>% group_by(Site, TimePt) %>% summarise(area = sum(Shape_Area)) %>%
  ggplot(aes(x = Site, y = area, fill = as.factor(TimePt))) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~TimePt)



#### Plot Google Site Map ####
  
setwd('C:/Users/Corinne.Amir/Documents/Github/ncrmp_common_maps')
utm = read.csv('data/misc/ncrmp_utm_zones.csv')                               # Load UTM zones (double check they're correct)
  

islands = c("GUA", "SAI", "TIN", "ROT", "MAU", "ASC", "FDP","AGR","AGU", "ALA","SAR","PAG","GUG") 
                       
utm_i = utm %>% subset(Island_Code %in%  islands)                                  # Subset required utm zone(s) (incorporate into function like Kisei later)
utm_i = utm %>% subset(Island_Code == islands)


load('data/gis_island_boundaries/ncrmp_islands_shp.RData')                    # Load island boundaries - shapefile vertices that have x,y coordinates
crs(ISL_bounds) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"    # Sets a coordinate reference system to the above shapefile --- not sure how/if this works


ISL_this = ISL_bounds[which(ISL_bounds$ISLAND_CD %in% islands),]  
  
ggmap::register_google("AIzaSyDpirvA5gB7bmbEbwB1Pk__6jiV4SXAEcY")

map <- get_map(location = c(mean(mari$Longitude, na.rm = T), mean(mari$Latitude, na.rm = T)), 
               maptype = "satellite",
               zoom = 6, 
               force = T)

mari_n <- mari %>% filter(Island == "MAU" | Island == "ASC")
map_n <- get_map(location = c(mean(mari_n$Longitude, na.rm = T), mean(mari_n$Latitude, na.rm = T)), 
                 maptype = "satellite",
                 zoom = 6, 
                 force = T)

mari_pag <- mari %>% filter(Island == "PAG" )
map_pag <- get_map(location = c(mean(mari_pag$Longitude, na.rm = T), mean(mari_pag$Latitude, na.rm = T)), 
                 maptype = "satellite",
                 zoom = 6, 
                 force = T)

mari_s <- mari %>% filter(Island == "GUA" | Island == "SAI")
map_s <- get_map(location = c(mean(mari_s$Longitude, na.rm = T), mean(mari_s$Latitude, na.rm = T)), 
                   maptype = "satellite",
                   zoom = 6, 
                   force = T)

dev.off()
mari_full <- ggmap(map) + 
             geom_polygon(data = ISL_this,
                          aes(long, lat, group=group),
                          fill = "darkgrey",
                          color = NA,
                          alpha = 0.9) +
             geom_label_repel(data = mari %>% distinct(Site, .keep_all = T), 
                              aes(Longitude, Latitude, label = Site),
                              box.padding = .8,
                              size = 4.5,
                              min.segment.length = .5,
                              segment.color = "white",
                              segment.size = 0.3,
                              nudge.x = .2) +
             geom_spatial_point(data = mari, aes(Longitude, Latitude), 
                                size = .5, 
                                color = "red", 
                                crs = 4326) +
             scale_x_continuous(limits=c(144.1,146.6)) +
             scale_y_continuous(limits=c(13.3,20.3)) 
             # annotation_scale(location = "br")
             # annotation_scale(location = "bl", width_hint = 0.3, height = unit(0.55, "cm"), text_cex = 1) + 

# ggsave("C:/Users/Corinne.Amir/Documents/Vital Rates/Analysis/MARAMP22/MARAMP22_map_full.pdf")
ggsave("C:/Users/Corinne.Amir/Documents/Vital Rates/Analysis/MARAMP22/MARAMP22_map.jpg",
       width = 5, height = 12)



mari_1 <- ggmap(map_n) + 
  geom_polygon(data = ISL_this,
               aes(long, lat, group=group),
               fill = "darkgrey",
               color = NA,
               alpha = 0.9) +
  geom_label_repel(data = mari_n %>% distinct(Site, .keep_all = T), 
                   aes(Longitude, Latitude, label = Site),
                   box.padding = .8,
                   size = 3.5,
                   min.segment.length = .7,
                   segment.color = "white",
                   segment.size = 0.3,
                   nudge.x = .2) +
  geom_spatial_point(data = mari_n, aes(Longitude, Latitude), 
                     size = 1.7, 
                     color = "red", 
                     crs = 4326) +
  scale_x_continuous(limits=c(145.1,145.5)) +
  scale_y_continuous(limits=c(19.6,20.1)) 
ggsave("C:/Users/Corinne.Amir/Documents/Vital Rates/Analysis/MARAMP22/MARAMP22_map_n.jpg")

mari_2 <- ggmap(map_pag) + 
  geom_polygon(data = ISL_this,
               aes(long, lat, group=group),
               fill = "darkgrey",
               color = NA,
               alpha = 0.9) +
  geom_label_repel(data = mari_pag %>% distinct(Site, .keep_all = T), 
                   aes(Longitude, Latitude, label = Site),
                   box.padding = .8,
                   size = 4,
                   min.segment.length = .7,
                   segment.color = "white",
                   segment.size = 0.3,
                   nudge.x = .2) +
  geom_spatial_point(data = mari_pag, aes(Longitude, Latitude), 
                     size = 3, 
                     color = "red", 
                     crs = 4326) +
  scale_x_continuous(limits=c(145.6,145.9)) +
  scale_y_continuous(limits=c(17.97,18.23)) 
ggsave("C:/Users/Corinne.Amir/Documents/Vital Rates/Analysis/MARAMP22/MARAMP22_map_pag.jpg")


mari_3 <- ggmap(map_s) + 
  geom_polygon(data = ISL_this,
               aes(long, lat, group=group),
               fill = "darkgrey",
               color = NA,
               alpha = 0.9) +
  geom_label_repel(data = mari_s %>% distinct(Site, .keep_all = T), 
                   aes(Longitude, Latitude, label = Site),
                   box.padding = 1.2,
                   size = 2.7,
                   min.segment.length = .8,
                   segment.color = "white",
                   segment.size = 0.3,
                   nudge.x = .2) +
  geom_spatial_point(data = mari_s, aes(Longitude, Latitude), 
                     size = .5, 
                     color = "red", 
                     crs = 4326) +
  scale_x_continuous(limits=c(144.4,145.8)) +
  scale_y_continuous(limits=c(13.2,15.4)) 



#### Format dataframe into archive csv file ####
# Add in Island_Code, DataorError, Error_Category
archive <- mari_col
colnames(mari)
head(archive)

archive <- archive %>% filter(Genus != "MOSP") # Remove MOSP because prevalence is too low


# Remove colonies <19cm2 in all time points
AdSize=(2.5^2*pi) 
`%notin%` <- Negate(`%in%`) 

t0 <- mari_col %>% filter(Year == "2014" & TL_Area < AdSize & Shape_Area < .0019) %>% 
  dplyr::select(c(Site, Genus, TL_Genet, Site_Genet))
t1 <- mari_col %>% filter(Year == "2017" & TL_Area < AdSize & Shape_Area < .0019) %>% 
  dplyr::select(c(Site, Genus, TL_Genet, Site_Genet)) 
t2 <- mari_col %>% filter(Year == "2022" & TL_Area < AdSize & Shape_Area < .0019) %>% 
  dplyr::select(c(Site, Genus, TL_Genet, Site_Genet))

step1 <- inner_join(t0,t1)
step1 <- step1 %>% distinct() 

step2 <- inner_join(t1,t2)
step2 <- step2 %>% distinct() 

smallcol <- rbind(step1,step2) %>% distinct()

archive<-subset(archive, Site_Genet %notin% smallcol$Site_Genet)


archive$Site <- sub("-", "_", archive$Site);archive$Site <- sub("-", "_", archive$Site) #twice for both underscores
archive$Error_Category <- "Growth Data"
archive$DataorError <- "DATA"
archive$Genet_full <- paste(archive$Site, archive$TL_Genet, archive$TimePt, sep = "_")
archive <- rename(archive, "Island_Code" = "Island")
archive <- archive %>% mutate(Island = case_when(Island_Code == "GUA" ~ "Guam",
                                                Island_Code == "MAU" ~ "Maug",
                                                Island_Code == "PAG" ~ "Pagan",
                                                Island_Code == "ASC" ~ "Asuncion",
                                                Island_Code == "SAI" ~ "Saipan"))
archive <- rename(archive, "Genus_Code" = "Genus")
archive <- archive %>% mutate(Genus = case_when(Genus_Code == "POSP" ~ "Porites sp.",
                                                Genus_Code == "MOSP" ~ "Montipora sp.",
                                                Genus_Code == "POCS" ~ "Pocillopora sp.",
                                                Genus_Code == "ACSP" ~ "Acropora sp."))
archive <- rename(archive, "Spec_Code" = "TL_Class")
archive <- rename(archive, "Date" = "TL_Date")
archive <- rename(archive, "ColonyName" = "Genet_full")
archive <- rename(archive, "Shape_Length" = "Shape_Leng")

archive <- archive %>% dplyr::select(-c(TL_Genet, TL_Cx, TL_Cy, TL_Area, TL_Perim, TL_Note, 
                                        Quadrat, Morph_Code, TimePt, Annotator, area_perim))
                    
archive <- archive %>% dplyr::select(Site, Island, Island_Code, Latitude, Longitude, Date, 
                              DataorError, Error_Category, ColonyName, Spec_Code, Genus,
                              Genus_Code, Shape_Length, Shape_Area)

head(archive)

#### Export Data ####

setwd('C:/Users/Corinne.Amir/Documents/GitHub/PIFSC_VitalRates/CSV files')
write.csv(mari,"ASRAMP23_VitalRates_patchlevel_CLEAN.csv",row.names = F)
write.csv(mari_col,"ASRAMP23_VitalRates_colonylevel_CLEAN.csv",row.names = F)
write.csv(archive,"ASRAMP23_VitalRates_colonylevel_archive.csv",row.names = F)

