library(tidyverse)

ll=load("./CSV files/Colony_Data_20210917_edited.rdata")
MA22=read.csv("./CSV files/ColonyTransitions/MARAMP22_ColonyTransitions.csv")
AS23=read.csv("./CSV files/ColonyTransitions/ASRAMP23_ColonyTransitions.csv")
names(AS23)
names(ColonyLevel)
HI19=ColonyLevel[,c("Island","Site","Genus_Code","ColonyID",
                    "StartingYear","StartingDate","StartingSize","StartingPerim",
                    "EndingYear","EndingDate","EndingSize","EndingPerim",
                    "TransitionTypeSimple","Interval_Years","log10_SS","log10_ES")]
HI19$log10_TransitionMagnitude=HI19$log10_ES-HI19$log10_SS
HI19$X=1:nrow(HI19)
names(MA22)
names(AS23)
names(HI19)
