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


#setwd('C:/Users/Corinne.Amir/Documents/GitHub/PIFSC_VitalRates/CSV files') # Github repo

raw <- read.csv("./Data/RawData/HARAMP24_VitalRates_12-02-2025.csv")
# raw <- read.csv("./CSV files/RawData/MARAMP22_VitalRates_06-24-2024.csv")
ll <- read.csv("./Data/MetaData/VitalRates_LatLong.csv")
effort <- read.csv("./Data/MetaData/VitalRates_SurveyEffort.csv")

#### QC Data ####
## (OPTIONAL) Remove superfluous columns:
colnames(raw)
#raw <- raw %>% select(-c(OID_, TL_SurfA))
raw <- raw %>% select(-TL_SurfA)

# Change some column names
raw <- raw %>% rename("Surface_Area" = "SArea")
                   

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
raw <- raw %>% mutate(Morph_Code = case_when(Morph_Code == "RM" ~ "EM",
                                             Morph_Code == "EM" ~ "EM",
                                             Morph_Code == "TB" ~ "TB",
                                             Morph_Code == "BR" ~ "BR",
                                             Morph_Code == "KN" ~ "KN",
                                             Morph_Code == "PL" ~ "PL",
                                             Morph_Code == "MD" ~ "MD"))
raw <- raw %>% mutate(Annotator = case_when (Annotator == "SJD" ~ "SD",
                                             Annotator == "SD" ~ "SD",
                                             Annotator == "CA" ~ "CA",
                                             Annotator == "ES" ~ "ES",
                                             Annotator == "KT" ~ "KT",
                                             Annotator == "MSL" ~ "ML"))


## Check if TimePt is labelled correctly:
a <- raw %>% group_by(Site, TL_Date, TimePt) %>% summarise(sum(TimePt)) ; View(a)
  # HOW-005-2017/2018 = Jan 1
  # TUT-019-2018 = Jan 1
# raw$TL_Date <- as.factor(raw$TL_Date)
# raw <- raw %>% mutate(TL_Date = recode(TL_Date,"1/1/2017" = "5/3/2017"))
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
vr %>% group_by(Site, TL_Date, Year) %>% summarise() # QC check



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
                                          TL_Class == "PLIC" ~ "POSP",    
                                          TL_Class == "PLUT" ~ "POSP",
                                          TL_Class == "POCS" ~ "POCS",
                                          TL_Class == "MOSP" ~ "MOSP",
                                          TL_Class == "MPAT" ~ "MOSP",
                                          TL_Class == "MCAP" ~ "MOSP",
                                          TL_Class == "ACSP" ~ "ACSP",
                                          TL_Class == "POSP" ~ "POSP"))
lapply(vr, unique) # double check


## Add island code
vr$Island <- str_sub(vr$Site,5,7)

## Add Lat and Long
ll <- ll %>% filter(Region != "vr") %>% rename(Site = ESD.Site.Name) %>% select(-Year)

vr <- left_join(vr, ll)


## Add m2 surveyed (collected from tracking spreadsheet)
vr$Year <- as.integer(vr$Year)
a <- left_join(vr, effort, by = c("Site", "Year","Genus"))


## Add area:perimeter ratio

vr$area_perim <- vr$Shape_Area/vr$Shape_Leng

head(vr)

## Add region code

# vr <- vr %>% mutate(Region = case_when(Island == "MAU" ~ "North",
#                                            Island == "PAG" ~ "North",
#                                            Island == "ASC" ~ "North",
#                                            Island == "SAI" ~ "South",
#                                            Island == "GUA" ~ "South"))

#### Consolidate into colony dataframe ####

# Create colony dataframe:
vr_col <- vr %>% dplyr::select(Genet_full, TL_Area, TL_Perim, Shape_Leng, Shape_Area, Surface_Area, area_perim) 
vr_meta <- vr %>% dplyr::select(Genet_full, Island, Site, TimePt, Year, TL_Date, Latitude,Longitude, #Region, Effort
                                Genus, TL_Class, TL_Genet, Quadrat) %>%
  distinct()
vr_col <- aggregate(.~Genet_full, data = vr_col, sum)
vr_col <- left_join(vr_meta,vr_col)

#group_by(Genet_full) %>% summarise(nrow(patches)) # Add in patch count 

# Total genets per Site-Year (ONLY colonies that survived 2+ time points):
vr_col$Site_Genet <- paste(vr_col$Site,vr_col$TL_Genet,  sep = "_") # use this column to filter


# Add column for number of patches
a <- vr %>% group_by(Site, Year,TL_Genet) %>% summarise(nPatches = n())

vr_col <- left_join(vr_col, a)


# Need to update/ consider minimum sizes....
#### Format dataframe into archive csv file #### 
# Add in Island_Code, DataorError, Error_Category
archive <- vr_col

# archive <- archive %>% filter(Genus != "MOSP") # Remove MARAMP22 MOSP because prevalence is too low


# Remove colonies <19cm2 in all time points
`%notin%` <- Negate(`%in%`) 

t0 <- vr_col %>% filter(Year == "2017" & Shape_Area < .0019) %>% 
  dplyr::select(c(Site, Genus, TL_Genet, Site_Genet))
t1 <- vr_col %>% filter(Year == "2018" & Shape_Area < .0019) %>% 
  dplyr::select(c(Site, Genus, TL_Genet, Site_Genet)) 
t2 <- vr_col %>% filter(Year == "2023" & Shape_Area < .0019) %>% 
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
                                                 Island_Code == "SAI" ~ "Saipan",
                                                 Island_Code == "BAK" ~ "Baker",
                                                 Island_Code == "HOW" ~ "Howland",
                                                 Island_Code == "OFU" ~ "Ofu & Olosega",
                                                 Island_Code == "ROS" ~ "Rose Atoll",
                                                 Island_Code == "TAU" ~ "Tau",
                                                 Island_Code == "TUT" ~ "Tutuila"))
archive <- rename(archive, "Genus_Code" = "Genus")
archive <- archive %>% mutate(Genus = case_when(Genus_Code == "POSP" ~ "Porites sp.",
                                                Genus_Code == "MOSP" ~ "Montipora sp.",
                                                Genus_Code == "POCS" ~ "Pocillopora sp.",
                                                Genus_Code == "ACSP" ~ "Acropora sp."))
archive <- archive %>% rename("Spec_Code" = "TL_Class",
                              "Date" = "TL_Date",
                              "ColonyName" = "Genet_full",
                              "Shape_Length" = "Shape_Leng")
archive <- archive %>% select(-c("TL_Area","TL_Perim"))

# archive <- archive %>% dplyr::select(-c(TL_Genet, TL_Area, TL_Perim,  
#                                         Quadrat, TimePt, area_perim))

archive <- archive %>% dplyr::select(Site, Island, Island_Code, Latitude, Longitude, Date, 
                                     DataorError, Error_Category, ColonyName, Spec_Code, Genus,
                                     Genus_Code, Shape_Length, Shape_Area, Surface_Area, nPatches)

head(archive)


#### Export Data ####

#setwd('C:/Users/Corinne.Amir/Documents/GitHub/PIFSC_VitalRates/CSV files')
write.csv(vr,"./Data/PatchLevel/HARAMP24_VitalRates_patchlevel_CLEAN.csv",row.names = F)
write.csv(vr_col,"./Data/ColonyLevel/HARAMP24_VitalRates_colonylevel_CLEAN.csv",row.names = F)
#write.csv(archive,"C:/Users/corinne.amir/Documents/Archiving/Vital Rates/ASRAMP/ASRAMP23_VitalRates_colonylevel_InportArchive.csv",row.names = F)
