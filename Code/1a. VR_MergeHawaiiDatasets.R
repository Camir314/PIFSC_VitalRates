# Incorporate older 2019 MHI and NWHI vital rates sites into current data structure
# Feb 2026

library(dplyr)
library(tidyr)
library(stringr)
library(anytime)


load("~/GitHub/PIFSC_VitalRates/Data/ColonyTransitions/HA_Trans/Colony_Data_20210917_edited.Rdata") # 2019 dataset
Col24=read.csv("./Data/ColonyTransitions/HARAMP24_ColonyTransitions.csv")
ll <- read.csv("./Data/MetaData/VitalRates_LatLong.csv")
effort <- read.csv("./Data/MetaData/VitalRates_SurveyEffort.csv")


# Remove sites resampled in 2024:
ColonyLevel %>% distinct(ColonyLevel$Site)
Col24 %>% distinct(Col24$Site)

rmsite = c("HAW_SIO_K08", "HAW_SIO_K10", "KUR_OCC_010", "LIS_OCC_005","MAI_SIO_K01","MAI_SIO_K02","FFS_OCC_014")

#`%notin%` <- Negate(`%in%`) 
ColonyLevel <- ColonyLevel %>% filter(Site %in% rmsite) 



# Update columns:
colnames(ColonyLevel)
colnames(Col24)

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
                                                       Site == "MAI_SIO_K02" ~ "SIO-MAI-K02")) %>%
                              select(-c("DataOrError","ColonyID","Spec_Code","SEC_NAME","REGION","SIG","SIS",
                                        "Genus","Fragmented","Recruit","Fragmentation","SiteInterval")) %>%
                              mutate(Interval = str_replace_all(Interval,"-","_")) %>%
                              mutate(Surface_Area_STA = NA_integer_,
                                     Surface_Area_END = NA_integer_) %>%
                              rename(Date_STA=StartingDate,Date_END=EndingDate,
                                     nPatches_STA=N_t0,nPatches_END=N_t1,
                                     Shape_Area_STA=StartingSize,Shape_Area_END=EndingSize,
                                     Shape_Leng_STA= StartingPerim,Shape_Leng_END=EndingPerim,
                                     Year_STA=StartingYear,Year_END=EndingYear)
                              
  
C <- ColonyLevel %>% select(Island,Site, Genus_Code, Site_Genet,Year_STA,Date_STA,Surface_Area_STA,
                            Shape_Area_STA,Shape_Leng_STA,nPatches_STA,Year_END,Date_END,Surface_Area_END,
                            Shape_Area_END,Shape_Leng_END,nPatches_END,TransitionTypeSimple,Interval)




## Add island code
vr$Island <- str_sub(vr$Site,5,7)

## Add Lat and Long
ll <- ll %>% filter(Region != "vr") %>% rename(Site = ESD.Site.Name) %>% select(-Year.Collected)

vr <- left_join(vr, ll)


## Add m2 surveyed
vr$Year <- as.integer(vr$Year)
a <- left_join(vr, effort, by = c("Site", "Year","Genus"))


## Add area:perimeter ratio

vr$area_perim <- vr$Shape_Area/vr$Shape_Leng


# # For converting the HA19 data to colony level rather than transition level data
# tdata <- ColonyLevel %>% select(c("uniqueID","Site","Island","Genus_Code","StartingDate","EndingDate",
#                                   "StartingPerim","EndingPerim","StartingSize","EndingSize", "N_t0", "N_t1")) #%>%
#                          # group_by(uniqueID) %>% 
#                          # pivot_longer(cols = StartingDate:EndingDate,
#                          #              names_to = "rm",
#                          #              values_to = "Date") %>%
#                          # select(-rm) %>%
#                          # pivot_longer(cols = c(StartingSize, EndingSize),
#                          #              names_to = "rm",
#                          #              values_to= "Shape_Area") %>%
#                          # select(-rm) 
# write.csv(tdata, "./Data/ColonyLevel/Archive/HARAMP19_VitalRates_patchlevel_export.csv",row.names = F) #reformatting in excel cause I give up
#  <- read.csv("./Data/Data/ColonyLevel/Archive/HARAMP19_VitalRates_patchlevel_import.csv")
# ColonyLevel2 = HARAMP19_VitalRates_patchlevel_import # re-import and make name more digestible 


# With re-imported haramp19 data, add in needed columns, create additional columns to match 2024 MHI/NWHI dataset

# ColonyLevel$Date <- anydate(ColonyLevel$Date) # Change date format

ColonyLevel$Year <- str_sub(ColonyLevel$Date,1,4) #create year column

ColonyLevel$Genet <- str_sub(ColonyLevel$uniqueID, -5,-1)

ColonyLevel$Site_Genet <- paste(ColonyLevel$Site,ColonyLevel$Genet,sep = "_")
                              
ColonyLevel2$Genet_full <- paste(ColonyLevel2$Site_Genet,ColonyLevel2$Year, sep = "_")

ColonyLevel2$area_perim <- ColonyLevel2$Shape_Area/ColonyLevel2$Shape_Leng

ll <- ll %>% rename("Site" = "ESD.Site.Name") %>% select(-Year.Collected) # add lat long values
ColonyLevel2 <- left_join(ColonyLevel2, ll)

ColonyLevel2 <- ColonyLevel2 %>% select(-uniqueID)

ColonyLevel2 <- rename(ColonyLevel2, "Genus" = "Genus_Code")                                        


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


