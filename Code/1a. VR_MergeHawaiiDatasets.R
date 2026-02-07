# Incorporate older 2019 MHI and NWHI vital rates sites into current data structure
# Feb 2026

library(dplyr)
library(tidyr)
library(stringr)



load("~/GitHub/PIFSC_VitalRates/Data/ColonyTransitions/HA_Trans/Colony_Data_20210917_edited.rdata") # 2019 dataset
#Colonylevel <- ColonyLevel19
load("~/GitHub/PIFSC_VitalRates/Data/ColonyTransitions/HA_Trans/Colony_Data_20210917_edited.rdata") # 2024 dataset

#HARAMP24_VitalRates_colonylevel_CLEAN <- ColonyLevel24


# Remove sites resampled in 2024:
ColonyLevel %>% distinct(ColonyLevel$Site)
HARAMP24_VitalRates_colonylevel_CLEAN %>% distinct(HARAMP24_VitalRates_colonylevel_CLEAN$Site)

rmsite = c("HAW_SIO_K08", "HAW_SIO_K10", "KUR_OCC_010", "LIS_OCC_005","MAI_SIO_K01","MAI_SIO_K02","FFS_OCC-002")

`%notin%` <- Negate(`%in%`) 
ColonyLevel <- ColonyLevel %>% filter(Site %in% rmsite) 


# Make site names match:
ColonyLevel <- ColonyLevel %>% mutate(Site = case_when(Site == "HAW_SIO_K08" ~ "SIO-HAW-K08",
                                                       Site == "HAW_SIO_K10" ~ "SIO-HAW-K10",
                                                       Site == "KUR_OCC_010" ~ "OCC-KUR-010",
                                                       Site == "LIS_OCC_005" ~ "OCC-LIS-005",
                                                       Site == "MAI_SIO_K01" ~ "SIO-MAI-K01",
                                                       Site == "MAI_SIO_K02" ~ "SIO-MAI-K02"))


# Update columns:
colnames(ColonyLevel)
colnames(HARAMP24_VitalRates_colonylevel_CLEAN)


# ColonyLevel$uniqueID <- paste(ColonyLevel$SIG,str_sub(ColonyLevel$ColonyID,20,24), sep = "_")
# a <- ColonyLevel %>% group_by(uniqueID) %>% filter(n()>1) %>% ungroup # double check unique identifier was achieved


tdata <- ColonyLevel %>% select(c("uniqueID","Site","Island","Genus","Genus_Code","StartingDate","EndingDate",
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
write.csv(tdata, "./Data/ColonyLevel/HARAMP19_VitalRates_patchlevel_export.csv",row.names = F) #reformatting in excel cause I give up



# With re-imported haramp19 data, add in needed columns

ColonyLevel2 = HARAMP19_VitalRates_patchlevel_import # make name more digestible 

ColonyLevel2$Year <- str_sub(ColonyLevel2$Date,-4,-1)

#create genet_full, timePt, latitude, longitude, area_perim, site_genet

                              
              
                                        

                              

