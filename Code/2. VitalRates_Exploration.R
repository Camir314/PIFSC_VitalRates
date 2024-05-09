#### MARAMP 2022 Vital Rates statistics and figure exploration
### Corinne Amir 
### Sept 2023

library(dplyr)
library(ggplot2)
library(ggmap)
library(viridis)
library(ggspatial)
library(ggrepel)
library(stringr)
library(sp)


# Make sure nfrags is followed through in the ColonyTransitions dataframe
# Match outputs from hawaii study IPM (colony names look similar)
# Can talk with just devin about how to take our data into her IPMs
# devin revisit and present code
# Get metadata and run through Kisei's EDS code (ask Juliette), link by sitevisit ID of each time we went
# Look into site target list for ASRAMP. Talk with Hannah/Carb budgets people to pick sites too. 
# Have Jonny/Frances start making those orthos and DEMs


# setwd("C:/Users/Corinne.Amir/Documents/N drive backup/Arc/Arc_Exports/")
setwd("C:/Users/Corinne.Amir/Documents/Vital Rates/Analysis/MARAMP22/CSV files/")
mari_col <- read.csv("MARAMP22_VitalRates_colonylevel_CLEAN.csv") # created in this script
# ll <- read.csv("MARAMP22_VitalRates_LatLong.csv")
# effort <- read.csv("MARAMP22_SurveyEffort.csv")
# mari <- read.csv("MARAMP22_VitalRates_patchlevel_CLEAN.csv") # created in this script
load("M:/FixedSiteProducts/Vital Rates Data/Patch_And_Colony_Data_20201103.rdata") # finished product from HARAMP19


#### Pre-treatments ####


  ##### Switch dataframe from long to wide ####
# Add in row for corals with recruitment or mortality

a <- mari_col %>% dplyr::select(Site_Genet, Year) %>% # Where does site_genet not occur in 2+ years?
  distinct() %>% 
  filter(Year != 2014) %>% # don't want to deal with 2014 yet
  add_count(Site_Genet) %>%
  filter(n>1)

aa <- anti_join(mari_col,a) 

dead <- aa %>% filter(Year == 2017) # Add rows for colonies that died by 2022
dead$Genet_full <- gsub('.{1}$' ,"1", dead$Genet_full)
dead$Year <-  rep(2022, times = length(dead$Year))
dead$TimePt <- rep(1, times = length(dead$TimePt))
dead$TL_Area <- rep(0, times = length(dead$TL_Area))
dead$Shape_Leng <- rep(0, times = length(dead$Shape_Leng))
dead$TL_Perim <- rep(0, times = length(dead$TL_Perim))
dead$Shape_Area <- rep(0, times = length(dead$Shape_Area))
dead$area_perim <- rep(0, times = length(dead$area_perim))

recruit <- aa %>% filter(Year == 2022) # Add rows for colonies that recruited in 2022
recruit$Genet_full <- gsub('.{1}$' ,"0", recruit$Genet_full)
recruit$Year <-  rep(2017, times = length(recruit$Year))
recruit$TimePt <- rep(0, times = length(recruit$TimePt))
recruit$TL_Area <- rep(0, times = length(recruit$TL_Area))
recruit$Shape_Leng <- rep(0, times = length(recruit$Shape_Leng))
recruit$TL_Perim <- rep(0, times = length(recruit$TL_Perim))
recruit$Shape_Area <- rep(0, times = length(recruit$Shape_Area))
recruit$area_perim <- rep(0, times = length(recruit$area_perim))

mari_col_rsd <- rbind(recruit, dead, mari_col) # combine recruit, survivor and dead data
mari_col_rsd <- mari_col_rsd %>% 
  dplyr::select(-c(Genet_full, TimePt, TL_Date,TL_Class, TL_Perim, area_perim, Shape_Leng))

t1 <- mari_col_rsd %>% filter(Year == 2017) %>% dplyr::select(-Year)
t2 <- mari_col_rsd %>% filter(Year == 2022) %>% dplyr::select(-Year)

t1 <- rename(t1, "Shape_Area_2017" = "Shape_Area",
                 "TL_Area_2017" = "TL_Area")
t2 <- rename(t2, "Shape_Area_2022" = "Shape_Area",
                 "TL_Area_2022" = "TL_Area")
t2$Shape_Area_2022 <- as.numeric(t2$Shape_Area_2022)
t1_t2 <- left_join(t1, t2)


# Deal with 2014 Maug sites separately 
a <- mari_col %>% filter(Island == "MAU") %>%
  dplyr::select(Site_Genet, Year) %>% # Where does site_genet not occur in each year?
  filter(Year != 2022) %>%
  distinct() %>% 
  add_count(Site_Genet) %>%
  filter(n!=2)
 
aa <- filter(mari_col,Site_Genet %in% a$Site_Genet) 

dead <- aa %>% filter(Year == 2014) # Add rows for colonies that died by 2017
dead$Genet_full <- gsub('.{1}$' ,"1", dead$Genet_full)
dead$Year <-  rep(2017, times = length(dead$Year))
dead$TimePt <- rep(1, times = length(dead$TimePt))
dead$TL_Area <- rep(0, times = length(dead$TL_Area))
dead$Shape_Leng <- rep(0, times = length(dead$Shape_Leng))
dead$TL_Perim <- rep(0, times = length(dead$TL_Perim))
dead$Shape_Area <- rep(0, times = length(dead$Shape_Area))
dead$area_perim <- rep(0, times = length(dead$area_perim))

dead2 <- dead # Make identical dataframe but for 2022 data
dead$Genet_full <- gsub('.{1}$' , "_2_", dead$Genet_full)
dead$Year <-  rep(2022, times = length(dead$Year))
dead$TimePt <- rep(2, times = length(dead$TimePt))

recruit <- aa %>% filter(Year == 2017) # Add rows for colonies that recruited in 2017
recruit$Genet_full <- gsub('.{1}$' ,"0", recruit$Genet_full)
recruit$Year <-  rep(2014, times = length(recruit$Year))
recruit$TimePt <- rep(0, times = length(recruit$TimePt))
recruit$TL_Area <- rep(0, times = length(recruit$TL_Area))
recruit$Shape_Leng <- rep(0, times = length(recruit$Shape_Leng))
recruit$TL_Perim <- rep(0, times = length(recruit$TL_Perim))
recruit$Shape_Area <- rep(0, times = length(recruit$Shape_Area))
recruit$area_perim <- rep(0, times = length(recruit$area_perim))

mau_col_rsd <- rbind(recruit, dead, dead2, mari_col) # combine recruit, survivor and dead data
mau_col_rsd <- mau_col_rsd %>% 
  dplyr::select(-c(Genet_full,TimePt, TL_Date,TL_Class, TL_Perim, area_perim, Shape_Leng))

t1 <- mau_col_rsd %>% filter(Year == 2014) %>% dplyr::select(-Year)
t2 <- mau_col_rsd %>% filter(Year == 2017) %>% dplyr::select(-Year)

t1 <- rename(t1, "Shape_Area_2014" = "Shape_Area",
                 "TL_Area_2014" = "TL_Area")
t2 <- rename(t2, "Shape_Area_2017" = "Shape_Area",
                 "TL_Area_2017" = "TL_Area")
t2$Shape_Area_2017 <- as.numeric(t2$Shape_Area_2017)
t1_t2_mau <- left_join(t1, t2)


col_wide <- full_join(t1_t2_mau, t1_t2) # Join 2014 dataframe with other dataframe
col_wide <- col_wide %>% select(Island, Site, Latitude, Longitude, Genus, TL_Genet,Quadrat, Effort, 
                                Site_Genet, Shape_Area_2014, TL_Area_2014, Shape_Area_2017, 
                                TL_Area_2017, Shape_Area_2022, TL_Area_2022)

col_wide$Shape_Area_2014 <- ifelse(col_wide$Island == "MAU" & is.na(col_wide$Shape_Area_2014), 0, as.numeric(col_wide$Shape_Area_2014))
col_wide$Shape_Area_2022 <- ifelse(is.na(col_wide$Shape_Area_2022), 0, as.numeric(col_wide$Shape_Area_2022))


# Add initial-final columns to wide data frame
# Add Interval_Years
a <- mari_col %>% group_by(Site, Year, TL_Date) %>%  
                  summarise(nColonies = n()) %>%
                  select(-nColonies)

b <- data.frame(Site = rep(unique(mari_col$Site), times = 3),
                   Year = rep(c("2014","2017","2022"), times = 10))

c <-left_join(b, a)                    

d <- spread(c, Year, TL_Date)
d <- d %>%  rename("Date_2014" = "2014",
                   "Date_2017" = "2017",
                   "Date_2022" = "2022") # Create separate columns for each time point

d$Interval_Years_1 <- as.numeric(difftime(d$Date_2017,d$Date_2014, units = "days")/365)
d$Interval_Years_2 <- as.numeric(difftime(d$Date_2022,d$Date_2017, units = "days")/365)

col_wide <- full_join(col_wide,d)

# Add TransitionMagnitude
col_wide$TransitionMagnitude_1 <- formatC((col_wide$Shape_Area_2017-col_wide$Shape_Area_2014),
                                         format = "f", digits = 6)
col_wide$TransitionMagnitude_2 <- formatC((col_wide$Shape_Area_2022-col_wide$Shape_Area_2017),
                                          format = "f", digits = 6)

col_wide$TransitionMagnitude_1 <- as.numeric(col_wide$TransitionMagnitude_1) 
col_wide$TransitionMagnitude_2 <- as.numeric(col_wide$TransitionMagnitude_2)


# Add TransitionRate
col_wide$TransitionRate_1 <- as.numeric(col_wide$TransitionMagnitude_1)/col_wide$Interval_Years_1
col_wide$TransitionRate_2 <- as.numeric(col_wide$TransitionMagnitude_2)/col_wide$Interval_Years_2


# Add TransitionType
  # Transition types used in HARAMP dataset: GROWTH,FUSION,SHRINK,FISSION, FUSION_FISSION MORT          
col_wide$TransitionType_1 <- ifelse(col_wide$Shape_Area_2014 == 0 & col_wide$Shape_Area_2017 > 0, "RECRUIT",
                             ifelse(col_wide$Shape_Area_2014 > 0 & col_wide$Shape_Area_2017 == 0, "MORT",
                             ifelse(col_wide$Shape_Area_2014 > 0 & col_wide$Shape_Area_2014 < col_wide$Shape_Area_2017, "GROWTH",
                             ifelse(col_wide$Shape_Area_2014 > 0 & col_wide$Shape_Area_2017 > 0 & col_wide$Shape_Area_2014 > col_wide$Shape_Area_2017, "SHRINK", NA))))
                                    
col_wide$TransitionType_2 <- ifelse(col_wide$Shape_Area_2017 == 0 & col_wide$Shape_Area_2022 > 0, "RECRUIT",
                             ifelse(col_wide$Shape_Area_2017 > 0 & col_wide$Shape_Area_2022 == 0, "MORT",
                             ifelse(col_wide$Shape_Area_2017 > 0 & col_wide$Shape_Area_2017 < col_wide$Shape_Area_2022, "GROWTH",
                             ifelse(col_wide$Shape_Area_2017 > 0 & col_wide$Shape_Area_2022 > 0 & col_wide$Shape_Area_2017 > col_wide$Shape_Area_2022, "SHRINK", NA))))     
                                  

# Add PercentChange
col_wide$PercentChange_1 <- ((col_wide$Shape_Area_2017-col_wide$Shape_Area_2014)/col_wide$Shape_Area_2014)*100
col_wide$PercentChange_2 <- ((col_wide$Shape_Area_2022-col_wide$Shape_Area_2017)/col_wide$Shape_Area_2017)*100


# Add Log2Ratio_Change
  # natural log


#### Figure Exploration ####
  #Scatter Plot 
  ggplot(t1_t2, aes(x = log(Shape_Area_2017), y = log(Shape_Area_2022), color = Region)) +
    geom_point(size = 1., alpha = .7) +
    geom_abline(size = 1) +
    scale_color_manual(values = alpha(c("blue", "red","green2","hotpink"),1)) + 
    scale_shape_manual(values = c(19)) + 
    facet_wrap(vars(Genus), nrow = 2, drop = F) +
    # facet_wrap(vars(Genus,Island), nrow = 3, drop = F) +
    theme_bw()
  
  ggplot(ColonyTransitions, aes(x = log(TL_Area_STA), y = log(TL_Area_END), color = Island)) +
    geom_point(size = 2., alpha = .7) +
    geom_abline(size = 1) +
    scale_color_manual(values = alpha(c("blue", "red","green2","hotpink","orange2"),1)) + 
    scale_shape_manual(values = c(19)) + 
    facet_wrap(vars(Genus), nrow = 2, drop = F) +
    # facet_wrap(vars(Genus,Island), nrow = 3, drop = F) +
    theme_bw()  
  

  # Box and whiskers plot: mean area
  ggplot(mari_col, aes(x = Region, y = log(Shape_Area), fill = as.factor(Year))) +
         geom_boxplot() +
         # geom_jitter(shape=16, position=position_jitter(0.2), color = as.factor(mari_col$Year)) +
         facet_wrap(vars(Genus))
  

  # Box and whiskers plot: perimeter to area ratio
  ggplot(mari_col, aes(x = Region, y = log(area_perim), fill = as.factor(Year))) +
    geom_boxplot() +
    facet_wrap(vars(Genus))
    # Looks the same as mean area box plot = double check area_perim calculation
    # Trying to show that colonies broke apart more in different regions/years 
  
  
  
  
  
  # Box and whiskers plot: proportional growth





#### Statistics Exploration ####

# Have to deal with spatial autocorrelation between sites on the same Islands
  # Maug sites are ~2km apart
  # Saipan sites are ~9km apart 
  # Guam sites are ~2/13/15km apart
  # Pagan sites are ~3km apart

# Have to deal with difference in sampling pressure at each site (survey package?)

# Environmental factors to look at: SST/DHW, Human population, Island type, Fish, Currents/island protection

  #### Kolmogorov-Smirnov tests ####
  # Break out individual sites, groups and years that look different based on kernel density plots/SFDs
  # Sites that have enough patches:
  # For POSP: all sites except GUA-015
  # For ACSP: maybe SAI-009
  # For POCS: MAU-019 (but 2017 only has three corals)
  
  # First, lets not split by site:
  posp0 <- mari %>%
    filter(Year %in% "2017",
           Genus == "POSP")
  posp1 <- mari %>%
    filter(Year %in% "2022",
           Genus == "POSP")
  
  acsp0 <- mari %>%
    filter(Year %in% "2017",
           Genus == "ACSP")
  acsp1 <- mari %>%
    filter(Year %in% "2022",
           Genus == "ACSP")
  
  pocs0 <- mari %>%
    filter(Year %in% "2017",
           Genus == "POCS")
  pocs1 <- mari %>%
    filter(Year %in% "2022",
           Genus == "POCS")
  
  # Second, split by populated vs unpopulated:
  # POSP
  unpop.posp0 <- mari %>% 
    filter(Island != "GUA" & Island != "SAI") %>%
    filter(Year %in% "2017",
           Genus == "POSP")
  unpop.posp1 <- mari %>% 
    filter(Island != "GUA" & Island != "SAI") %>%
    filter(Year %in% "2022",
           Genus == "POSP")
  
  pop.posp0 <- mari %>% 
    filter(Island == "GUA" | Island == "SAI") %>%
    filter(Year %in% "2017",
           Genus == "POSP")
  pop.posp1 <- mari %>% 
    filter(Island == "GUA" | Island == "SAI") %>%
    filter(Year %in% "2022",
           Genus == "POSP")
  
  #ACSP
  unpop.acsp0 <- mari %>% 
    filter(Island != "GUA" & Island != "SAI")%>%
    filter(Year %in% "2017",
           Genus == "ACSP") 
  unpop.acsp1 <- mari %>% 
    filter(Island != "GUA" & Island != "SAI") %>%
    filter(Year %in% "2022",
           Genus == "ACSP") 
  
  pop.acsp0 <- mari %>% 
    filter(Island == "GUA" | Island == "SAI") %>%
    filter(Year %in% "2017",
           Genus == "ACSP") # Guam only contributes 1 coral
  pop.acsp1 <- mari %>% 
    filter(Island == "GUA" | Island == "SAI") %>%
    filter(Year %in% "2022",
           Genus == "ACSP") # Guam only contributes 1 coral
  
  # POCS
  unpop.pocs0 <- mari %>% 
    filter(Island != "GUA" & Island != "SAI")%>%
    filter(Year %in% "2017",
           Genus == "POCS") # only 27 corals, none from Pagan
  unpop.pocs1 <- mari %>% 
    filter(Island != "GUA" & Island != "SAI") %>%
    filter(Year %in% "2022",
           Genus == "POCS") # Pagan only contributes 2 corals
  
  pop.pocs0 <- mari %>% 
    filter(Island == "GUA" | Island == "SAI") %>%
    filter(Year %in% "2017",
           Genus == "POCS") 
  pop.pocs1 <- mari %>% 
    filter(Island == "GUA" | Island == "SAI") %>%
    filter(Year %in% "2022",
           Genus == "POCS") # only 26 corals
  
  # Kolmogorov-Smirnov tests:
  ks.test(posp0$TL_Area, posp1$TL_Area) # D = 0.064877, p-value = 0.009682
  ks.test(acsp0$TL_Area, acsp1$TL_Area) # D = 0.22441, p-value = 0.0002861
  ks.test(pocs0$TL_Area, pocs1$TL_Area) # D = 0.1688, p-value = 0.04478
  
  ks.test(unpop.posp0$TL_Area, unpop.posp1$TL_Area) # D = 0.094853, p-value = 0.004505
  ks.test(unpop.acsp0$TL_Area, unpop.acsp1$TL_Area) # D = 0.38083, p-value = 0.0001103
  ks.test(unpop.pocs0$TL_Area, unpop.pocs1$TL_Area) # D = 0.56427, p-value = 8.995e-07
  
  ks.test(pop.posp0$TL_Area, pop.posp1$TL_Area) # D = 0.206, p-value = 2.83e-11
  ks.test(pop.acsp0$TL_Area, pop.acsp1$TL_Area) # D = 0.31988, p-value = 0.001423
  ks.test(pop.pocs0$TL_Area, pop.pocs1$TL_Area) # D = 0.32456, p-value = 0.01164
  
  ks.test(pop.posp0$TL_Area, unpop.posp0$TL_Area) # D = 0.23539, p-value < 2.2e-16
  ks.test(pop.acsp0$TL_Area, unpop.acsp0$TL_Area) # D = 0.18831, p-value = 0.1432: ACSP distribution NOT different between populated and unpopulated islands before bleaching event!
  ks.test(pop.pocs0$TL_Area, unpop.pocs0$TL_Area) # D = 0.50975, p-value = 6.377e-05
  
  ks.test(pop.posp1$TL_Area, unpop.posp1$TL_Area) # D = 0.10806, p-value = 0.002973
  ks.test(pop.acsp1$TL_Area, unpop.acsp1$TL_Area) # D = 0.58855, p-value = 3.7e-10
  ks.test(pop.pocs1$TL_Area, unpop.pocs1$TL_Area) # D = 0.37908, p-value = 0.0004607
  
  
  #### Anovas ####
  # Create dataframes
  posp <- mari %>% filter(Genus == "POSP")
  acsp <- mari %>% filter(Genus == "ACSP")
  pocs <- mari %>% filter(Genus == "POCS")
  
  
  
  # Is mean size different between Islands 
  # Assumptions
  bartlett.test(TL_Area~as.factor(Year), data = posp) # significant: unequal variances
  bartlett.test(TL_Area~as.factor(Year), data = acsp) # significant
  bartlett.test(TL_Area~as.factor(Year), data = pocs) # significant
  tapply(posp$TL_Area, list(as.factor(posp$Year)), shapiro.test) # all significant: non-normal distribution
  
  bartlett.test(TL_Area~interaction(as.factor(OBS_YEAR), DEPTH_BIN), data = swa) # significant
  tapply(swa$AVG_DEPTH, list(swa$OBS_YEAR), shapiro.test) # all significant
  tapply(swa$AVG_DEPTH, list(swa$DEPTH_BIN), shapiro.test) # all significant
  
  
  #two-way anova
  summary(aov(CCA~as.factor(Year)*TRT, data = cca)) #trt and year are significant but interaction is not
  TukeyHSD(aov(CCA~as.factor(Year), data = cca)) 
  # 1709 > 1510, 1606, 1609
  TukeyHSD(aov(CCA~TRT, data = cca)) 
  # REM > BYSP, RAN, XSP
  