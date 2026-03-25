# Incorporate older 2019 MHI and NWHI vital rates sites into current data structure
# Feb 2026

library(dplyr)
library(tidyr)
library(stringr)
library(anytime)


load("~/GitHub/PIFSC_VitalRates/Data/ColonyTransitions/HA_Trans/Colony_Data_20210917_edited.Rdata") # 2019 dataset
CT24=read.csv("./Data/ColonyTransitions/HARAMP24_ColonyTransitions.csv")
CL24=read.csv("./Data/ColonyLevel/HARAMP24_VitalRates_colonylevel_CLEAN.csv")
ll <- read.csv("./Data/MetaData/VitalRates_LatLong.csv")
effort <- read.csv("./Data/MetaData/VitalRates_SurveyEffort.csv")


# Remove sites resampled in 2024:
ColonyLevel %>% distinct(ColonyLevel$Site)
CL24 %>% distinct(CL24$Site)

rmsite = c("HAW_SIO_K08", "HAW_SIO_K10", "KUR_OCC_010", "LIS_OCC_005","MAI_SIO_K01","MAI_SIO_K02","FFS_OCC_014")

`%notin%` <- Negate(`%in%`)
ColonyLevel <- ColonyLevel %>% filter(Site %in% rmsite) 



# Update columns:
colnames(ColonyLevel)
colnames(CL24)

ColonyLevel$Site_Genet <- paste(ColonyLevel$Site,str_sub(ColonyLevel$ColonyID,20,24), sep = "_")
ColonyLevel$uniqueID <- paste(ColonyLevel$SIG,str_sub(ColonyLevel$ColonyID,20,24), sep = "_") #create unique ID#
ColonyLevel %>% group_by(uniqueID) %>% filter(n()>1) %>% ungroup # double check unique identifier was achieved


# Make site names match:
ColonyLevel <- ColonyLevel %>% mutate(Site = case_when(Site == "FFS_OCC_014" ~ "OCC-FFS-014",
                                                       Site == "HAW_SIO_K08" ~ "SIO-HAW-K08",
                                                       Site == "HAW_SIO_K10" ~ "SIO-HAW-K10",
                                                       Site == "KUR_OCC_010" ~ "OCC-KUR-010",
                                                       Site == "LIS_OCC_005" ~ "OCC-LIS-005",
                                                       Site == "MAI_SIO_K01" ~ "SIO-MAI-K01",
                                                       Site == "MAI_SIO_K02" ~ "SIO-MAI-K02")) 
ColonyLevel$Island <- str_sub(ColonyLevel$Site,5,7)

dup_col <- ColonyLevel %>% group_by(Site_Genet) %>% # QC check 
  filter(n() > 1) %>%
  ungroup()
head(dup_col) 
  
# duplicate_counts <- ColonyLevel  %>%
#   group_by_all() %>%
#   summarise(count = n()) %>%
#   filter(count > 1)

#######################
# For converting HA19 into colony transitions (skipping colony level) - still run [partially] through script 2
#######################
ColonyLevel <- ColonyLevel %>% mutate(Interval = str_replace_all(Interval,"-","_")) %>%
                               mutate(Surface_Area_STA = NA_integer_,
                                      Surface_Area_END = NA_integer_) %>%
                               rename(Date_STA=StartingDate,Date_END=EndingDate,
                                      nPatches_STA=N_t0,nPatches_END=N_t1,
                                      Shape_Area_STA=StartingSize,Shape_Area_END=EndingSize,
                                      Shape_Leng_STA= StartingPerim,Shape_Leng_END=EndingPerim,
                                      Year_STA=StartingYear,Year_END=EndingYear)
                              
  
Col19 <- ColonyLevel %>% select(Island, Site, Genus_Code, Site_Genet,Year_STA,Date_STA,Surface_Area_STA,
                            Shape_Area_STA, Shape_Leng_STA, nPatches_STA, Year_END, Date_END, Surface_Area_END,
                            Shape_Area_END, Shape_Leng_END, nPatches_END, TransitionTypeSimple, Interval,
                            Interval_Years, Mortality, Survival);head(C)


write.csv(Col19,"./Data/ColonyTransitions/HARAMP19_VitalRates_ColonyTransitions_NotUpdated.csv",row.names = F)


#########################
# For converting the HA19 data to colony level rather than transition level data
# Finished reformatting colony transition -> colony level in excel cause I give up
########################

tdata <- ColonyLevel %>% select(c("uniqueID","Site_Genet","Site","Island","Genus_Code","StartingDate","EndingDate",
                                  "StartingPerim","EndingPerim","StartingSize","EndingSize", "N_t0", "N_t1")) #%>%
                         # group_by(uniqueID) %>%
                         # pivot_longer(cols = StartingDate:EndingDate,
                         #              names_to = "rm",
                         #              values_to = "Date") %>%
                         # select(-rm) %>%
                         # pivot_longer(cols = c(StartingSize, EndingSize),
                         #              names_to = "rm",
                         #              values_to= "Shape_Area") %>%
                         # select(-rm)
write.csv(tdata, "./Data/ColonyLevel/Archive/HARAMP19_VitalRates_colonylevel_export.csv",row.names = F) # file I need to edit

CL19 <- read.csv("./Data/ColonyLevel/Archive/HARAMP19_VitalRates_colonylevel_import.csv") # file I edited = use this to finish QC
# CL19$Island <- str_sub(CL19$Site,5,7) # error missed when importing back in

# Found colonies where size_end of transition 1 doesn't match size_start of transition 2 for the same colony ID:
rmsitegenet = c("FFS_OCC_014_00530","FFS_OCC_014_00534","HAW_SIO_K10_01880","LIS_OCC_005_02142","MAI_SIO_K01_02363",
                "MAI_SIO_K01_02396","MAI_SIO_K01_02433","MAI_SIO_K01_02439","MAI_SIO_K01_02450","MAI_SIO_K01_02482",
                "MAI_SIO_K02_02648","MAI_SIO_K02_02810","MAI_SIO_K02_02845","MAI_SIO_K02_02901","MAI_SIO_K02_02981",
                "MAI_SIO_K02_03001") 
CL19 <- CL19 %>% filter(Site_Genet %notin% rmsitegenet) %>% # remove these colonies and duplicates from the dataset
                 select(-uniqueID) %>%
                 distinct()


##########################
## With re-imported haramp19 data, add in columns, create additional columns to match 2024 MHI/NWHI colony-level dataset
##########################

ll <- ll %>% select(-c(Island,Year.Collected))
CL19 <- left_join(CL19, ll, by = "Site") # Add Lat and Long

CL19$Date <- anydate(CL19$Date) # Change date format

CL19$Year <- str_sub(CL19$Date,1,4) #create year column
CL19$MoYear <- paste(str_sub(CL19$Date,6,7),CL19$Year, sep = "-")

CL19$Genet <- str_sub(CL19$Site_Genet, -5,-1)
                              
CL19$Genet_full <- paste(CL19$Site_Genet,CL19$Year, sep = "_") # unique ID

CL19$area_perim <- CL19$Shape_Area/CL19$Shape_Leng # Add area:perimeter ratio

#CL19 <- CL19 %>% select(-uniqueID) %>% distinct()

head(CL19)
dup_col <- CL19 %>% group_by(Genet_full) %>% # QC check 
  filter(n() > 1) %>%
  ungroup()
unique(dup_col$Site_Genet)


## Add m2 surveyed
ColonyLevel$Year <- as.integer(ColonyLevel$Year)
a <- left_join(ColonyLevel, effort, by = c("Site", "Year","Genus"))

# ## Add m2 surveyed UPDATE!
# vr$Year <- as.integer(vr$Year)
# a <- left_join(vr, effort, by = c("Site", "Year","Genus"))
# # Add annotator
# # Add morph code (not in original csv file)


# Remove dummy rows added for mortality/recruitment
## Look for rows with NA values
rows_with_na <- ColonyLevel2 %>%
  filter(if_any(everything(), is.na))
print(rows_with_na) # only includes rows where length and area = 0 = ok but remove these rows to match 2024 dataset

ColonyLevel2 <- ColonyLevel2 %>% filter(Shape_Area != 0)



# Repeating same QC from script 1

## Make sure Genus matches for all colonies with the same ID#
diff_species <- ColonyLevel2 %>% group_by(Genet_full, Genus) %>% 
  summarise() %>% 
  filter(n() >1) %>% 
  ungroup()
print(diff_species)

diff_species_diffyear <- ColonyLevel2 %>% group_by(Site_Genet, Genus) %>% 
  summarise() %>% 
  filter(n() >1) %>% 
  ungroup()
print(diff_species_diffyear) # for now, only going to edit where genus isn't equal


# Export data 
write.csv(ColonyLevel2,"./Data/ColonyLevel/HARAMP19_VitalRates_colonylevel_CLEAN.csv",row.names = F)


