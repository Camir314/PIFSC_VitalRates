rm(list=ls())
# Libraries ---------------------------------------------------------------
library(igraph)
library(lubridate)
library(ggplot2)
library(ggpubr)
library(plyr)
library(stringr)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(flextable)
library(officer)


# ## Functions ------------------------------------------------------------

#Pass dataset of single site, multiple timepoint data.frame "Patches"
#Must Have first column as "PatchID"
Patches2ColonyGraphs=function(Patches){
  require(igraph)
  require(lubridate)
  # Error Handling -------------------------------------------------
  #Force PatchName to first column
  #patch first
  if(which(names(Patches)=="PatchName")!=1){stop("Patches2ColonyGraphs: ERROR: 'Patches' data frame must have first column with unique Patch Names  named 'PatchName'...")}
  #need cols
  neededcols=c("LinkTP_For","LinkTP_Bac","Date","TimePt","Circrat","PatchID","PatchName","Notes","Shape_Area")#,"FID"
  missingcols=setdiff(neededcols,names(Patches))
  if(length(missingcols)>0){
    missingcols=setdiff(neededcols,names(Patches))
    stop(paste0("Patches2ColonyGraphs: ERROR: 'Patches' data frame must have columns: ",paste0(missingcols,collapse=",")))
  }
  # How Many Digits in PatchID -------------------------------------------------
  NIDdigit=mean(nchar(Patches$PatchID)-2) #assuming "P_" then leading zero padded number for patchID
  if(any(nchar(Patches$PatchID)!=(NIDdigit+2))){
    stop(paste0("Patches2ColonyGraphs: ERROR: This code assumes that PatchIDs start with 'P_' and then have a leading-zero-padded number. You've got variable length IDs, which is going to lead to issues."))
  }
  #CleanDates
  Patches$R_DATE=mdy(Patches$Date)
  #Add Type
  Patches$Type="LIVE"
  
  # Build To/From DataFrame -------------------------------------------------
  Transitions=NULL
  print(paste("Building Transitions from",nrow(Patches),"patches."))
  for (i in 1:nrow(Patches)){
    #PatchID-->PatchName
    P_Home=Patches$PatchName[i]
    PN_sp=unlist(strsplit(x=P_Home,split="_P_"))
    PN_pre=paste0(PN_sp[1],"_")
    PN_post=paste0("_",unlist(strsplit(PN_sp[2],"_"))[2])
    P_To=unlist(strsplit(x = as.character(Patches$LinkTP_For[i]),split = ","))
    P_From=unlist(strsplit(x = as.character(Patches$LinkTP_Bac[i]),split = ","))
    if(!all(is.na(P_To))){#If there's anything aside from NA in "P_To"
      P_To=na.omit(P_To)
      P_To[which(P_To!="-99")]=paste0(PN_pre,"P_",leadz(P_To,n=NIDdigit),PN_post)
      P_To[which(P_To=="-99")]=paste0("MORT_",P_Home)
      T_To=data.frame(From=P_Home,To=P_To,stringsAsFactors = F)
      Transitions=rbind(Transitions,T_To)
    }
    if(!all(is.na(P_From))){#If there's anything aside from NA in "P_From"
      P_From=na.omit(P_From)
      P_From[which(P_From!="-99")]=paste0(PN_pre,"P_",leadz(P_From,n=NIDdigit),PN_post)
      P_From[which(P_From=="-99")]=paste0("RECR_",P_Home)
      T_From=data.frame(From=P_From,To=P_Home,stringsAsFactors = F)
      Transitions=na.omit(rbind(Transitions,T_From))
    }
  }
  #Shoudl expect many repeated to-from links
  Transitions=unique(Transitions)
  print(paste("Generated",nrow(Transitions),"unique transitions from",nrow(Patches),"patches."))
  
  #   #Add Patches for -99s (recruitment/mortality) -------------------------
  #Figure out patches to add (names)
  AllPatches=sort(unique(c(Transitions[,1],Transitions[,2])))
  N99raw=setdiff(AllPatches,Patches$PatchName)
  
  #Here NPP should flag any transition that is either RECR MORT or some breakdown between PatchID and PatchName
  #This last group can include (1) Trans species links, and so other issues. For now, I'm deleting the non MORT/RECR transitions, reporting and moving on 
  Problem_i=setdiff(1:length(N99raw),c(grep(pattern = "MORT",x=N99raw),grep(pattern = "RECR",x=N99raw)))
  
  if(length(Problem_i)>0){
    print(paste("Dropping",length(Problem_i),"Transitions with problematic PatchNames from",nrow(Transitions),"total."))
    print(paste("Examples include:",paste(N99raw[sample(Problem_i,20,replace = T)],collapse=", ")))
    
    #Set N99
    N99=N99raw[-Problem_i]
    #Remove from Transitions
    Problem_Trow=sort(unique(c(which(Transitions[,1]%in%N99raw[Problem_i]),which(Transitions[,2]%in%N99raw[Problem_i]))))
    if(length(Problem_Trow)>0){Transitions=Transitions[-Problem_Trow,]}
  }else{
      N99=N99raw
  }
                
  if(length(N99)>0){
    N99ID=substr(x = N99,start = 6,stop = nchar(N99))
    N99TYPE=substr(x = N99,start = 1,stop = 4)
    
    #Build and modify DF
    N99DF=Patches[match(N99ID,Patches$PatchName),]
    #Get TimePoint LookUp
    TPlu=unique(N99DF[,c("Site","TimePt","Date","R_DATE")])
    TPlu$STP=paste0(TPlu$Site,"_",TPlu$TimePt)
    
    #Modify Patch Data - PatchID
    N99DF$PatchName=N99
    #Modify Patch Data - PatchID
    N99DF$PatchID=substr(N99ID,start = 18,nchar(N99ID)-2)
    #Modify Patch Data - FID
    N99DF$FID=max(Patches$FID,na.rm=T)+1:nrow(N99DF)
    #Modify Patch Data - Site
    N99DF$Site=substr(N99ID,1,11)
    #Modify Patch Data - Circat
    N99DF$Circrat=NA
    #Modify Patch Data - Shape_Area
    N99DF$Shape_Area=0
    #Modify Patch Data - LinkTP_Bac
    N99DF[which(N99TYPE=="MORT"),"LinkTP_Bac"]=N99ID[which(N99TYPE=="MORT")]
    N99DF[which(N99TYPE=="RECR"),"LinkTP_Bac"]=NA
    #Modify Patch Data - LinkTP_For
    N99DF[which(N99TYPE=="RECR"),"LinkTP_For"]=N99ID[which(N99TYPE=="RECR")]
    N99DF[which(N99TYPE=="MORT"),"LinkTP_For"]=NA
    #Modify Patch Data - TimePt
    N99DF[which(N99TYPE=="MORT"),"TimePt"]=N99DF[which(N99TYPE=="MORT"),"TimePt"]+1
    N99DF[which(N99TYPE=="RECR"),"TimePt"]=N99DF[which(N99TYPE=="RECR"),"TimePt"]-1
    #Modify Patch Data - Date
    N99STP=paste0(N99DF$Site,"_",N99DF$TimePt)
    N99DF$Date=TPlu[match(N99STP,TPlu$STP),"Date"]
    N99DF$R_DATE=TPlu[match(N99STP,TPlu$STP),"R_DATE"]
    #Modify Patch Data - Notes
    N99DF$Notes=N99TYPE
    N99DF$Type=N99TYPE
    
    Patches_rm=rbind(Patches,N99DF)
  }else{
    Patches_rm=Patches
  }
  
  g1 <- graph_from_data_frame(Transitions, directed=TRUE, vertices=Patches_rm)
  
  E(g1)$AreaChange_cm2=10^4*(head_of(g1,E(g1))$Shape_Area-tail_of(g1,E(g1))$Shape_Area)
  E(g1)$AreaChangeSign=sign(E(g1)$AreaChange_cm2)
  
  #Set Transition Types
  E(g1)$TRANSITION_TYPE=NA
  E(g1)$TRANSITION_TYPE[E(g1)$AreaChangeSign>=0]="GROWTH"
  E(g1)$TRANSITION_TYPE[E(g1)$AreaChangeSign<0]="SHRINK"
  
  MORT_i=which(head_of(g1,E(g1))$Type=="MORT")
  E(g1)$TRANSITION_TYPE[MORT_i]="MORT"
  RECR_i=which(tail_of(g1,E(g1))$Type=="RECR")
  E(g1)$TRANSITION_TYPE[RECR_i]="RECR"
  
  #FUsion/FIssion are slightly harder because they apply in Edge-index space, but assess in Vertex-index space
  #Mark Every Edge From a vertex with >1 out going connections a "FISSION" edge
  FISSION_Vi=which(degree(g1,v=unique(tail_of(graph = g1,es=E(g1))),
                          mode="out")>1)#index in V() space
  E(g1)[.from(names(FISSION_Vi))]$TRANSITION_TYPE="FISSION"
  #Mark Every Edge From a vertex with >1 in-coming connections a "FUSION" edge
  FUSION_Vi=which(degree(g1,
                         v=unique(head_of(graph = g1,es=E(g1))),
                         mode="in")>1)
  E(g1)[.to(names(FUSION_Vi))]$TRANSITION_TYPE="FUSION"
  
  #Mark Every Edge Listed as both FUSION and FISSION a FUSION_FISSION Edge
  FUSFIS_Ei=intersect(E(g1)[.from(names(FISSION_Vi))], E(g1)[.to(names(FUSION_Vi))])
  E(g1)[FUSFIS_Ei]$TRANSITION_TYPE="FUSION_FISSION"
  
  E(g1)$Interval_Years=(head_of(g1,E(g1))$R_DATE-tail_of(g1,E(g1))$R_DATE)/365.25
  E(g1)$DataOrError=as.vector(cut(E(g1)$Interval_Years,breaks=c(0,7/365.25,90/365.25,999),labels=c("ERROR","NEITHER","DATA"),include.lowest=T))
  E(g1)$CrossAnnotator=head_of(g1,E(g1))$Annotation_Round!=tail_of(g1,E(g1))$Annotation_Round
  E(g1)$GrowthRate_cm2_yr=E(g1)$AreaChange_cm2/E(g1)$Interval_Years
  
  #Assign ColonyID to Vertex
  V(g1)$ColonyID="NA"
  sg1=decompose.graph(g1,min.vertices = 2)
  for(icol in 1:length(sg1)){
    thisgraph=sg1[[icol]]
    V(g1)[which(V(g1)$name%in%V(thisgraph)$name)]$ColonyID=paste0("C_",formatC(icol,width=5,flag=0),"_",V(thisgraph)$Site[1])
  }
  return(g1)
}

#Plotting Support
NodeTypeColor=c("green","red","pink")
names(NodeTypeColor)=c("LIVE","MORT","RECR")

EdgeTypeColor=c("darkblue","lightblue","brown","yellow","red","pink","orange")
names(EdgeTypeColor)=c("GROWTH","SHRINK","FUSION","FISSION","MORT","RECR","FUSION_FISSION")

SizeScaleFunc=function(x,lowest=5,highest=15,trans="none"){
  if(trans=="log"){
    x_l=log(x+1)
    x_sc=(x_l-min(x_l,na.rm=T))/(max(x_l,na.rm=T)-min(x_l,na.rm=T))
  }else{
    x_sc=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}
  sz=lowest+x_sc*(highest-lowest)
  return(sz)
}

subg_Vname=function(graph,vertex_names){
  sg1 <- decompose.graph(graph,mode="weak")
  subfollow=function(s){if(any(V(s)$name %in% vertex_names)) V(s)$name else NULL}  
  neighverts <- unique(unlist(sapply(sg1,FUN=subfollow)))
  gout <- induced.subgraph(graph=graph,vids=neighverts)
  return(gout)
}

subg_Vattr=function(graph,attr,comparison,INFO=T){
  sg1 <- decompose.graph(graph,mode="weak")
  eval(parse(text=paste0("vertex_names=V(graph)$name[which(V(graph)$",attr,comparison,")]")))
  subfollow=function(s){if(any(V(s)$name %in% vertex_names)) V(s)$name else NULL}  
  neighverts <- unique(unlist(sapply(sg1,FUN=subfollow)))
  gout <- induced.subgraph(graph=graph,vids=neighverts)
  if(INFO){print(paste0(length(decompose.graph(gout,mode="weak"))," subgraphs of ",length(sg1)," with feature identified."))}
  return(gout)
}

subg_Eind=function(graph,edge_indices){
  verts=V(graph)[.inc(E(graph)[edge_indices])]$name
  return(subg_Vname(graph,verts))
}

as.nv=function(x){return(as.numeric(as.vector(x)))}
leadz=function(x,n){return(formatC(as.numeric(as.vector(x)),width=n,flag=0))} #formats numbers individually
#Converts value to vector and then converts to numeric vector
#if the value is a charcter/word, returns NAs and gives warning message NAs introduced by coercion

mod1=function(x,k){return(1+mod(x-1,k))}
plot.colgraph=function(graph,NodeColor=NodeTypeColor,EdgeColor=EdgeTypeColor,SizeFunc=SizeScaleFunc,layoutfun=layout_nicely){
  shapes=c("circle", "square","pie","sphere", "rectangle","crectangle", "vrectangle")
  #eval(parse(text=paste0("glay=",layoutfun,"(graph)")))
  glay=layoutfun(graph)
  plot.igraph(graph,#layout=glay,
              edge.width=SizeFunc(E(graph)$AreaChange_cm2,lowest=.05,highest=6),
              edge.arrow.width=SizeFunc(abs(E(graph)$AreaChange_cm2),lowest=.05,highest=4),
              edge.arrow.size=SizeFunc(abs(E(graph)$AreaChange_cm2),lowest=.05,highest=1),
              edge.color=EdgeColor[E(graph)$TRANSITION_TYPE],
              vertex.color=NodeColor[V(graph)$Type],
              vertex.pie.color=NodeColor[V(graph)$Type],
              vertex.size=SizeFunc(V(graph)$Shape_Area,lowest=3,highest=10,trans="log"),
              vertex.shape=shapes[mod1(V(graph)$TimePt,length(shapes))],
              vertex.label.cex=0.75,
              vertex.label.color="black",
              vertex.label.dist=1.5)
}

# Load Datasets -----------------------------------------------------------
basepath="/Users/c-rod/Documents/GitHub/ArcMap_Exports/"  
#basepath="M:/FixedSiteProducts/Vital Rates Data/ArcMap Exports-20200903T211716Z-001/ArcMap Exports/"
#"C:/Users/Thomas.Oliver/WORK/Projects/Vital Rates/PatchToTransitionCode/" #Tom's dataset
#basepath="/Users/c-rod/Documents/GitHub/ArcMapExports_ICRS/" #Tom's dataset



#Set up to loop through files
fl=list.files(path = basepath,pattern = ".csv",full.names = T)
#Need to distinguish annotator by file
annotation_num=rep(1,length(fl))
annotation_num[grep(pattern = "VitalRatesQC",x =fl)]=2

CanonColNames = c("FID","Site","Date","TimePt","Annotation_Round","Circrat","Empty_Skip_Circrat_Flag","PatchID","PatchName",
                  "Spec_Code","Morph_Code","BLE_EXT","BLE_SEV","Shape_Length","Shape_Area","Shape_Diameter","LinkTP_Bac","LinkTP_For","Notes","FileName")
# extra link columns need becasue they ran out of room, need to merege

#ColonyID is the Empty flag..
#Site _ PatchID should be unique ##### Dups are likely due to double annotatotion

#load 'em all in storing different columnnames/orders in differentlists
for(i_f in 1:length(fl)){
  thissite=read.csv(fl[i_f],stringsAsFactors = FALSE)
  thissite$Annotation_Round=annotation_num[i_f]
  thissite$FileName=strsplit(fl[i_f],"/")[[1]][6]
  
  names(thissite)
#Catch easy naming errors
  if(any(names(thissite)%in%"ColonyID")){names(thissite)[which(names(thissite)=="ColonyID")]="Empty_Skip_Circrat_Flag"}
  if(any(names(thissite)%in%"OBJECTID")){names(thissite)[which(names(thissite)=="OBJECTID")]="FID"}
  if(any(names(thissite)%in%"SpeciesCod")){names(thissite)[which(names(thissite)=="SpeciesCod")]="Spec_Code"}
  if(any(names(thissite)%in%"Poly_Area")){names(thissite)[which(names(thissite)=="Poly_Area")]="Shape_Area"}
  if(any(names(thissite)%in%"Poly_Perim")){names(thissite)[which(names(thissite)=="Poly_Perim")]="Shape_Length"}
  if(any(names(thissite)%in%"Shape_Leng")){names(thissite)[which(names(thissite)=="Shape_Leng")]="Shape_Length"}
  if(any(names(thissite)%in%"MBG_Diamet")){names(thissite)[which(names(thissite)=="MBG_Diamet")]="Shape_Diameter"}
  if(any(names(thissite)%in%"Max_Diam")){names(thissite)[which(names(thissite)=="Max_Diam")]="Shape_Diameter"}
  if(any(names(thissite)%in%"Max_Diameter")){names(thissite)[which(names(thissite)=="Max_Diameter")]="Shape_Diameter"}
  names(thissite)
  
  #Take care of NA in EMPTY
  thissite$Empty_Skip_Circrat_Flag[is.na(thissite$Empty_Skip_Circrat_Flag)]="NOT_EMPTY"
  thissite$Empty_Skip_Circrat_Flag[str_trim(thissite$Empty_Skip_Circrat_Flag)==""]="NOT_EMPTY"
  
  #Drop any BS extra columns
  dropcols=c("Id","Shape_Area.1","Shape_Leng.1","MBG_Diamet","ORIG_FID","Area","MBG_Diameter")
  thissite=thissite[,setdiff(names(thissite),dropcols)]
  names(thissite)
  
  #If there's no bleaching or MorphCode or Diameter data, add NA cols
  if(!"BLE_EXT"%in%names(thissite)){thissite$BLE_EXT=NA}
  if(!"BLE_SEV"%in%names(thissite)){thissite$BLE_SEV=NA}
  if(!"Morph_Code"%in%names(thissite)){thissite$Morph_Code=NA}
  if(!"Shape_Diameter"%in%names(thissite)){thissite$Shape_Diameter=NA}
  names(thissite)
  
  ###Tidy Link Columns:
  #Get blank and NA cells to null string ""
  LinkCols_i=grep(pattern = "LinkTP_",names(thissite))
  for(lc_i in 1:length(LinkCols_i)){thissite[which(thissite[,LinkCols_i[lc_i]]==" "|is.na(thissite[,LinkCols_i[lc_i]])),LinkCols_i[lc_i]]=""}
  
  #Concatenate multiple link columns
  if(length(grep(pattern = "LinkTP_F",names(thissite)))>1){
    print(paste(i_f,":",fl[i_f],names(thissite)[grep(pattern = "LinkTP_F",names(thissite))]))
    LF=paste(thissite$LinkTP_For,thissite$LinkTP_F_1,sep=",")
    LF[which(LF==",")]=""
    thissite$LinkTP_For=LF
    thissite=thissite[,setdiff(names(thissite),"LinkTP_F_1")]
  }
  if(length(grep(pattern = "LinkTP_B",names(thissite)))>1){
    print(paste(i_f,":",fl[i_f],names(thissite)[grep(pattern = "LinkTP_B",names(thissite))]))
    LB=paste(thissite$LinkTP_Bac,thissite$LinkTP_B_1,sep=",")
    LB[which(LB==",")]=""
    thissite$LinkTP_Bac=LB
    thissite=thissite[,setdiff(names(thissite),"LinkTP_B_1")]
  }
  #Make sure all periods are commas, drop timepoint after -99s
  perstr_iF=grep(thissite$LinkTP_For,pattern="\\.")
  perstr_iB=grep(thissite$LinkTP_Bac,pattern="\\.")
  #for each period in a string, either turn any -99.X into just -99s, or change any periods into commas
  #For For
  if(length(perstr_iF)>0){
    for(i_ps in 1:length(perstr_iF)){
      cellstr=thissite$LinkTP_For[perstr_iF[i_ps]]
      cellstr=gsub(cellstr,pattern="\\..",replacement="\\.")
      thisps=strsplit(cellstr,",")
      tpsnum=lapply(thisps,as.numeric)
      #from list find/replace in strings the -99s
      l1_i=which(lapply(thisps,length)==1)
      if(length(l1_i)>0) {neg99s_i=l1_i[tpsnum[[l1_i]]<0]}else{neg99s_i=NULL}
      if(length(neg99s_i)>0){thisps[[neg99s_i]]="-99"} else {print(paste0("Problem with this csv file: ", fl[i_f], " at line ", perstr_iF[i_ps], " for these problematic linkages: ", cellstr))}
      #moosh it back together, then replace any . with ,
      thisstr=paste(unlist(thisps),collapse=",")
      thissite$LinkTP_For[perstr_iF[i_ps]]=gsub(thisstr,pattern="\\.",replacement="\\,")
    }
  }
  
  #For Bac
  if(length(perstr_iB)>0){
    for(i_ps in 1:length(perstr_iB)){
      cellstr=thissite$LinkTP_Bac[perstr_iB[i_ps]]
      cellstr=gsub(cellstr,pattern="\\..",replacement="\\.")
      thisps=strsplit(cellstr,",")
      tpsnum=lapply(thisps,as.numeric)
      #from list find/replace in strings the -99s
      l1_i=which(lapply(thisps,length)==1)
      if(length(l1_i)>0) {neg99s_i=l1_i[tpsnum[[l1_i]]<0]}else{neg99s_i=NULL}
      if(length(neg99s_i)>0){thisps[[neg99s_i]]="-99"} else {print(paste0("Problem with this csv file: ", fl[i_f], " at line ", perstr_iB[i_ps], " for these problematic linkages: ", cellstr))}
      #moosh it back together, then replace any . with ,
      thisstr=paste(unlist(thisps),collapse=",")
      thissite$LinkTP_Bac[perstr_iB[i_ps]]=gsub(thisstr,pattern="\\.",replacement="\\,")
    }
  }
  
  #Prep for final addition
  thissite=thissite[,sort(names(thissite))]
  
  #Here's where we actually add the new data into our datastructure
  if(i_f==1){ # if this is the first time we're doing this, just set up a list for the data frame and one for the column names
    formatlist=list(names(thissite))
    stacklist=list(thissite)
  }else{ #if this isnt the first time, it gets a bit more complicated...
    formatmatch=NULL
    for(i_for in 1:length(formatlist)){#for each element in the format list
      #check if the current dataframe's columnames match the list-stored sets of columnames
      formatmatch=c(formatmatch,all(formatlist[[i_for]]==names(thissite)))
    }
    formi=which(formatmatch)
    if(length(formi)==0){#There is no match format existing, New Foramt!
      print(paste0("New Format! i_f =",i_f,". Total of ",length(formatlist)+1))
      formatlist[[length(formatlist)+1]]=names(thissite)
      stacklist[[length(stacklist)+1]]=thissite
    }else{
      stacklist[[formi]]=rbind(stacklist[[formi]],thissite)
    }
  }

  print(paste0(i_f," of ",length(fl)))
  #write out THISSITE WITH dATE CODE
  
}
length(stacklist) #want stacklist to be 1 b/c this means all the formats match
lapply(stacklist,nrow)

formatlist
Patches=stacklist[[1]][,CanonColNames]

#exclude the pseudo circrats that are used to calculate area surveyed
dim(Patches)
Patches <- Patches %>% filter(Empty_Skip_Circrat_Flag == "NOT_EMPTY")
#there are 4 patches w/ no site, date, patch link, etc. Remove them
Patches <- Patches %>% filter(PatchID>0)
dim(Patches)

Patches$R_DATE=mdy(Patches$Date)
#write out PaTCHES WITH dATE CODE

#Clean Up And Re-Org PatchID
Pnames=names(Patches) 
Pnames=Pnames[-which(Pnames=="PatchID")]
Patches$PatchIDraw=Patches$PatchID 
Patches$PatchID=paste0("P_",leadz(Patches$PatchIDraw,n=5))#,"_",Patches$Site,"_",Patches$Annotation_Round)
Patches=Patches[,c("PatchID",Pnames)] 

uS=             c("MAI_SIO_K01", "FFS_OCC_002","FFS_OCC_014","HAW-OCC-002","HAW_OCC_003",
                  "HAW_OCC_010","HAW_SIO_K08","HAW_SIO_K10","KUR_OCC_010",
                  "OCC_MAI_00","MAI_SIO_K02","MAI_SIO_OL3","OAH-OCC-010","OAH22",
                  "PHR-OCC-016","OAH_D2","D2_8in","MAI-002","OAH-22","OAH23","OAH_SIO_99","SpiralTest","SpilTestA_8in","SpiTest_AB_6A")
#sites need to be in alphabetical order
SiteCoherenceLU=c("FFS_OCC_002","FFS_OCC_014","HAW_OCC_002","HAW_OCC_003",
                  "HAW_OCC_010","HAW_SIO_K08","HAW_SIO_K10","KUR_OCC_010",
                  "LIS_OCC_005","MAI_OCC_002","MAI_SIO_K01","MAI_SIO_K02",
                  "MAI_SIO_OL3","OAH_OCC_010","OAH_XXX_022",
                  "PHR_OCC_016") 
#removed: "OAH_D2x_xxx","OAH_D2x_xxx","OAH_SPI_0AB","OAH_SPI_00A","OAH_SPI_00A","OAH_SPI_0AB",
names(SiteCoherenceLU)=unique(Patches$Site) 
Patches$Site=SiteCoherenceLU[Patches$Site] # if SiteCoherence isn't in alphabetical order all of the Sites get shuffled to the wrong PatchName b/c they share the same PatchID

table(Patches$Site)

Patches$PatchName=paste(Patches$Site,Patches$Spec_Code,Patches$PatchID,Patches$Annotation_Round,sep="_")
#Clean Up And Re-Org PatchName
Pnames=names(Patches)
Pnames=Pnames[-which(Pnames=="PatchName")]
Patches=Patches[,c("PatchName",Pnames)]
PatchesWithEmpty=Patches
Patches=subset(PatchesWithEmpty,Empty_Skip_Circrat_Flag=="NOT_EMPTY")

GenusCodeLUnames=c("PLIG","PMEA","PGRA","POCS","POSP","PLUT","PLIC","PLOB","MFLA","MPAT","MOSP","MCAP")
GenusCodeLU=     c("POCS","POCS","POCS","POCS","POSP","POSP","POSP","POSP","MOSP","MOSP","MOSP","MOSP")
names(GenusCodeLU)=GenusCodeLUnames

Patches$Spec_Code=str_trim(Patches$Spec_Code)
drop_SpBlank=which(Patches$Spec_Code=="")
if(length(drop_SpBlank)>0){Patches=Patches[-drop_SpBlank,]}
Patches$Genus_Code=GenusCodeLU[Patches$Spec_Code]
table(Patches$Genus_Code,useNA="always")

BadDates_i=c(which(year(Patches$R_DATE)<2000),which(year(Patches$R_DATE)>2020))
table(Patches$Site[BadDates_i])


# Patches2ColonyGraphs Call -----------------------------------------------

#For every Ann_Round==2 Patch, match it to the linked patches from Ann_Round==1
A2_i=which(Patches$Annotation_Round==2)
#Check Annotator's Patch Notes
AndNote=sort(str_trim(unique(Patches$Notes[A2_i])))
GoodNotes=AndNote[c(1,2,5,6,9,22,23,25,26,27)]
BadNotes=setdiff(AndNote,GoodNotes)
BadNote_i=which(Patches$Notes%in%BadNotes)
BadNote_PatchName=Patches$PatchName[BadNote_i]

#MatchUp Patch Up
PN2=Patches$PatchName[A2_i]
PN=substr(PN2,1,nchar(Patches$PatchName[A2_i])-2)
PN1=paste0(PN,"_1")

#Match 'em up
A1_im=match(PN1,Patches$PatchName)

#Check for Link Columns that don't match between 1 and 2
#  BacMM1=A1_im[which(Patches$LinkTP_Bac[A2_i]!=Patches$LinkTP_Bac[A1_im])]
#  ForMM1=A1_im[which(Patches$LinkTP_For[A2_i]!=Patches$LinkTP_For[A1_im])]
#  BacMM2=A2_i[which(Patches$LinkTP_Bac[A2_i]!=Patches$LinkTP_Bac[A1_im])]
#  ForMM2=A2_i[which(Patches$LinkTP_For[A2_i]!=Patches$LinkTP_For[A1_im])]
#  cbind(Patches[BacMM1,c("PatchName","LinkTP_Bac")],Patches[BacMM2,c("PatchName","LinkTP_Bac")])
#  cbind(Patches[ForMM1,c("PatchName","LinkTP_For")],Patches[ForMM2,c("PatchName","LinkTP_For")])
# Patches$LinkTP_Bac[BacMM2]=Patches$LinkTP_Bac[BacMM1]
# Patches$LinkTP_For[ForMM2]=Patches$LinkTP_For[ForMM1]
# table(Patches$LinkTP_For[A2_i]==Patches$LinkTP_For[A1_im])
# table(Patches$LinkTP_Bac[A2_i]==Patches$LinkTP_Bac[A1_im])
# #There are lots, because Annotator only deliniated "Error" time points once in most sites, ...
# table(Patches$Annotation_Round,Patches$TimePt,Patches$Site)
#just 1: HAW_OCC_002, HAW_OCC_003,HAW_OCC_010,HAW_SIO_K08,HAW_SIO_K10,KUR_OCC_010,
# MAI_OCC_002

#Patches to Colony

CGall=Patches2ColonyGraphs(Patches = Patches)

VitalRate_Growth=
    data.frame(Site=head_of(CGall,E(CGall))$Site,
               DataOrError=E(CGall)$DataOrError,
               Annotator_Tail=tail_of(CGall,E(CGall))$Annotation_Round,
               Annotator_Head=head_of(CGall,E(CGall))$Annotation_Round,
               ColonyID=head_of(CGall,E(CGall))$ColonyID,
               T0_PatchName=tail_of(CGall,E(CGall))$name,
               T1_PatchName=head_of(CGall,E(CGall))$name,
               EdgeLabels=paste0(tail_of(CGall,E(CGall))$name,"-",head_of(CGall,E(CGall))$name),
               Spec_Code=tail_of(CGall,E(CGall))$Spec_Code,
               Genus_Code=tail_of(CGall,E(CGall))$Genus_Code,
               AreaSurveyed=0.5*length(unique(V(CGall)$Circrat)),
               StartingDate=as.Date(tail_of(CGall,E(CGall))$R_DATE,origin="1970/01/01"),
               EndingDate=as.Date(head_of(CGall,E(CGall))$R_DATE,origin="1970/01/01"),
               Interval_Years=E(CGall)$Interval_Years,
               StartingSize=10^4*tail_of(CGall,E(CGall))$Shape_Area,
               EndingSize=10^4*head_of(CGall,E(CGall))$Shape_Area,
               StartingPerimeter=10^2*tail_of(CGall,E(CGall))$Shape_Length,
               EndingPerimeter=10^2*head_of(CGall,E(CGall))$Shape_Length,
               StartingMaxDiam=10^2*tail_of(CGall,E(CGall))$Shape_Diameter,
               EndingMaxDiam=10^2*head_of(CGall,E(CGall))$Shape_Diameter,
               TransitionMagnitude=E(CGall)$AreaChange_cm2,
               TransitionRate=E(CGall)$GrowthRate_cm2_yr,
               TransitionType=E(CGall)$TRANSITION_TYPE)
VitalRate_Growth$PercentChange=100*VitalRate_Growth$TransitionMagnitude/VitalRate_Growth$StartingSize
VitalRate_Growth$Log2Ratio_Change=log2(VitalRate_Growth$EndingSize/VitalRate_Growth$StartingSize)
head(VitalRate_Growth)


save(file = "/Users/c-rod/Documents/GitHub/Patch-To-Transition/ThesisSites_2021_09_CR.rdata",list=c("Patches","CGall","VitalRate_Growth"))


# Load and plot -----------------------------------------------------------


load(file = "/Users/c-rod/Documents/GitHub/Patch-To-Transition/ThesisSites_2021_07_CR.rdata")
GenusCodeLUnames=c("PMEA","PGRA","POSP","PLUT", "PLIC","MFLA","MPAT", "MOSP","POCS","POSP","PLOB","MCAP")
GenusCodeLU=c("POCS","POCS","POSP","POSP","POSP", "MOSP","MOSP", "MOSP","POCS","POSP","POSP","MOSP")
names(GenusCodeLU)=GenusCodeLUnames


#Convert Species Based PatchName to Genus Based Patch Name-T0
VitalRate_Growth$T0_PatchName=as.character(VitalRate_Growth$T0_PatchName)
VitalRate_Growth$T0_PatchName_G=paste0(substr(VitalRate_Growth$T0_PatchName,1,12),
                                       GenusCodeLU[substr(VitalRate_Growth$T0_PatchName,13,16)],
                                       substr(VitalRate_Growth$T0_PatchName,17,nchar(VitalRate_Growth$T0_PatchName)))
RM_i=which(substr(VitalRate_Growth$T0_PatchName,1,4)%in%c("RECR","MORT"))
VitalRate_Growth$T0_PatchName_G[RM_i]=paste0(substr(VitalRate_Growth$T0_PatchName[RM_i],1,17),
                                       GenusCodeLU[substr(VitalRate_Growth$T0_PatchName[RM_i],18,21)],
                                       substr(VitalRate_Growth$T0_PatchName[RM_i],22,nchar(VitalRate_Growth$T0_PatchName[RM_i])))

#Convert Species Based PatchName to Genus Based Patch Name-T1
VitalRate_Growth$T1_PatchName=as.character(VitalRate_Growth$T1_PatchName)
VitalRate_Growth$T1_PatchName_G=paste0(substr(VitalRate_Growth$T1_PatchName,1,12),
                                       GenusCodeLU[substr(VitalRate_Growth$T1_PatchName,13,16)],
                                       substr(VitalRate_Growth$T1_PatchName,17,nchar(VitalRate_Growth$T1_PatchName)))
RM_i=which(substr(VitalRate_Growth$T1_PatchName,1,4)%in%c("RECR","MORT"))
VitalRate_Growth$T1_PatchName_G[RM_i]=paste0(substr(VitalRate_Growth$T1_PatchName[RM_i],1,17),
                                             GenusCodeLU[substr(VitalRate_Growth$T1_PatchName[RM_i],18,21)],
                                             substr(VitalRate_Growth$T1_PatchName[RM_i],22,nchar(VitalRate_Growth$T1_PatchName[RM_i])))

VRG_GS=subset(VitalRate_Growth,TransitionType%in%c("GROWTH","SHRINK"))
VRG_GSerr=subset(VRG_GS,DataOrError=="ERROR")


#Cross_Annotator Error: DATA, ERROR, CROSS_SAMETP, CROSS_DIFFTP
#Patches$PatchName=paste(Patches$Site,Patches$Spec_Code,Patches$PatchID,Patches$Annotation_Round,sep="_")
#Genus-Based Patch Name
PN_G=paste(Patches$Site,Patches$Genus_Code,Patches$PatchID,Patches$Annotation_Round,sep="_")
S2=Patches$Spec_Code[which(Patches$Annotation_Round==2)]
G2=Patches$Genus_Code[which(Patches$Annotation_Round==2)]

#Identify Patches from Annotator 2, build equivalent Genus PatchName for Annotator 1
crossPNXS=data.frame(Spec_Code=S2,Genus_Code=G2,PN_0=PN_G[which(Patches$Annotation_Round==2)],stringsAsFactors = F)
crossPNXS$PN_1=paste0(substr(crossPNXS$PN_0,1,nchar(crossPNXS$PN_0)-2),"_1")
crossPNXS$DataOrError="DIFF.AN_SAME.IM"

#Drop No-Match Patches
drop=which(!crossPNXS$PN_1%in%PN_G);
print(paste("Dropping",length(drop),"non-existent patches from",nrow(crossPNXS),"candidates."));
if(length(drop)>0){crossPNXS=crossPNXS[-drop,]}

#Given Ann2 Patches, Find Ann1 Patches, Match Ann2 with Linked: Track TO-1 to linked T1-1 patch, and T1-1 to T0-1
crossPNXX=data.frame(Spec_Code=crossPNXS$Spec_Code,Genus_Code=crossPNXS$Genus_Code,PN_0=crossPNXS$PN_0,stringsAsFactors = F)
PN_T0erri=which(crossPNXS$PN_1%in%VRG_GSerr$T0_PatchName_G)
PN_T1erri=which(crossPNXS$PN_1%in%VRG_GSerr$T1_PatchName_G)
crossPNXX$PN_1=NA
crossPNXX$PN_1[PN_T0erri]=VRG_GSerr$T1_PatchName_G[match(crossPNXS$PN_1[PN_T0erri],VRG_GSerr$T0_PatchName_G)]
crossPNXX$PN_1[PN_T1erri]=VRG_GSerr$T0_PatchName_G[match(crossPNXS$PN_1[PN_T1erri],VRG_GSerr$T1_PatchName_G)]
crossPNXX$DataOrError="DIFF.AN_DIFF.IM"

#Drop mismatches
drop=which(!crossPNXX$PN_1%in%PN_G);
print(paste("Dropping",length(drop),"non-existent patches from",nrow(crossPNXX),"candidates."));
if(length(drop)>0){crossPNXX=crossPNXX[-drop,]}

#Join 'em up
crossPN=rbind(crossPNXS,crossPNXX)
drop=which(!crossPN$PN_0%in%PN_G);print(paste("Dropping",length(drop),"non-existent patches"));if(length(drop)>0){crossPN=crossPN[-drop,]}
drop=which(!crossPN$PN_1%in%PN_G);print(paste("Dropping",length(drop),"non-existent patches"));if(length(drop)>0){crossPN=crossPN[-drop,]}

#Merge CrossPN with VRG_GSerr
#names(VitalRate_Growth)
ErrorFile.=VRG_GS[,c("Site","DataOrError","T0_PatchName_G","T1_PatchName_G","Spec_Code","Genus_Code","StartingSize","EndingSize","TransitionMagnitude")]
addCross=data.frame(Site=substr(crossPN$PN_0,1,11),
                    DataOrError=crossPN$DataOrError,
                    T0_PatchName_G=crossPN$PN_0,
                    T1_PatchName_G=crossPN$PN_1,
                    Spec_Code=crossPN$Spec_Code,
                    Genus_Code=crossPN$Genus_Code,
                    StartingSize=10^4*Patches$Shape_Area[match(crossPN$PN_0,PN_G)],
                    EndingSize=10^4*Patches$Shape_Area[match(crossPN$PN_1,PN_G)])
addCross$TransitionMagnitude=addCross$EndingSize-addCross$StartingSize
ErrorFile=rbind(ErrorFile.,addCross)

#Randomize Annotator 1 and 2
RAND_ANN=FALSE
if(RAND_ANN){
  err_i=which(ErrorFile$DataOrError!="DATA")
  swap_i=err_i[which(sample(x = c(T,F),size = length(err_i),replace = T))]
  tmpPNt0=ErrorFile$T0_PatchName_G[swap_i]
  ErrorFile$T0_PatchName_G[swap_i]=ErrorFile$T1_PatchName_G[swap_i]
  ErrorFile$T1_PatchName_G[swap_i]=tmpPNt0
  tmpStSize=ErrorFile$StartingSize[swap_i]
  ErrorFile$StartingSize[swap_i]=ErrorFile$EndingSize[swap_i]
  ErrorFile$EndingSize[swap_i]=tmpStSize
}

#Calculate Delta Metrics
ErrorFile$TransitionMagnitude=(ErrorFile$EndingSize)-(ErrorFile$StartingSize)
ErrorFile$PercentChange=100*(ErrorFile$EndingSize-ErrorFile$StartingSize)/ErrorFile$StartingSize
ErrorFile$Log2Ratio=log2(ErrorFile$EndingSize/ErrorFile$StartingSize)
ErrorFile$LogNRatio=log(ErrorFile$EndingSize/ErrorFile$StartingSize)
ErrorFile$Log10Ratio=log10(ErrorFile$EndingSize/ErrorFile$StartingSize)
ErrorFile$DeltaLog2Area=log2(ErrorFile$EndingSize)-log2(ErrorFile$StartingSize)
ErrorFile$DeltaLogNArea=log(ErrorFile$EndingSize)-logN(ErrorFile$StartingSize)
ErrorFile$DeltaLog10Area=log10(ErrorFile$EndingSize)-log10(ErrorFile$StartingSize)
ErrorFile=subset(ErrorFile,Genus_Code%in%c("POCS","POSP"))

#Clean UP DataOrError
ErrorFile$DataOrError=as.character(as.vector(ErrorFile$DataOrError))
ErrorFile$DataOrError[ErrorFile$DataOrError=="ERROR"]="SAME.AN_DIFF.IM"
ErrorFile$Error_Category=as.vector(ErrorFile$DataOrError)
ErrorFile$Error_Category[ErrorFile$Error_Category%in%c("DATA")]="Growth Data"
ErrorFile$Error_Category[ErrorFile$Error_Category%in%c("DIFF.AN_SAME.IM")]="Cross-Annotator Error"
ErrorFile$Error_Category[ErrorFile$Error_Category%in%c("DIFF.AN_DIFF.IM")]="Imaging & Annotator Error"
ErrorFile$Error_Category[ErrorFile$Error_Category%in%c("SAME.AN_DIFF.IM")]="Imaging Error"
ErrorFile$Error_Category=factor(ErrorFile$Error_Category,levels=c("Growth Data","Imaging Error","Cross-Annotator Error","Imaging & Annotator Error"))

#Clean up Genus
ErrorFile$Genus=as.vector(ErrorFile$Genus_Code)
ErrorFile$Genus[ErrorFile$Genus=="MOSP"]="Montipora sp."
ErrorFile$Genus[ErrorFile$Genus=="POCS"]="Pocillopora sp."
ErrorFile$Genus[ErrorFile$Genus=="POSP"]="Porites sp."
ErrorFile$Genus=factor(ErrorFile$Genus,levels=c("Pocillopora sp.","Porites sp.","Montipora sp."))

#Flag Bad Annotator Notes
BadNote_PatchName_G=paste0(substr(BadNote_PatchName,1,12),GenusCodeLU[substr(BadNote_PatchName,13,16)],substr(BadNote_PatchName,17,nchar(BadNote_PatchName)))
ErrorFile$AnnotatorNoteFlag=FALSE
ErrorFile$AnnotatorNoteFlag[which(ErrorFile$T0_PatchName_G%in%BadNote_PatchName_G)]=TRUE

#Restrict to Adults (5 cm diam)
AdSize=(2.5^2*pi)
AdErrorFile=subset(ErrorFile,StartingSize>AdSize&EndingSize>AdSize)

#select Metric 
AdErrorFile$Y=AdErrorFile$TransitionMagnitude
YPlotLabel="Colony Size Change (cm^2)"
outname="AbsChange"

#Describe Bounds (with outliers)
RMSE=function(x){mn=mean(x,na.rm=T)
return(sqrt(mean((mn-x)^2,na.rm=T)))}

oSC=3
adbound=ddply(AdErrorFile,.(Genus,Error_Category),summarize,
              N=length(Y),
              MN=mean(Y),
              MD=median(Y),
              SD=sd(Y),
              SE=sd(Y)/sqrt(length(Y)),
              CI95=1.96*SE,
              RMSE=RMSE(Y),
              LOW3SD=MN-oSC*SD,HIGH3SD=MN+oSC*SD,
              LOW95=MN-CI95,HIGH95=MN+CI95,
              TMBias_p=t.test(Y,mu=0)$p.val)

#Define outlier
AdErrorFile$Outlier=NA
for(i in 1:nrow(adbound)){
  thisset_i=which(AdErrorFile$Genus==adbound$Genus[i]&AdErrorFile$Error_Category==adbound$Error_Category[i])
  AdErrorFile$Outlier[thisset_i]=(AdErrorFile$Y[thisset_i]<adbound$LOW3SD[i]|AdErrorFile$Y[thisset_i]>adbound$HIGH3SD[i])
}
DROP_BADNOTE=F
if(DROP_BADNOTE){
  AdErrorFile$Outlier[which(AdErrorFile$AnnotatorNoteFlag)]=TRUE
}
table(AdErrorFile$Outlier,AdErrorFile$DataOrError)

#Describe Bounds (without outliers)
AdErrorFileNO=subset(AdErrorFile,Outlier==F)
AEFYES=subset(AdErrorFile,Outlier==T)
adboundNO=ddply(AdErrorFileNO,.(Genus,Error_Category),summarize,
                N=length(Y),
                MN=mean(Y),
                MD=median(Y),
                SD=sd(Y),
                SE=sd(Y)/sqrt(N),
                CI95=1.96*SE,
                RMSE=RMSE(Y),
                LOW3SD=MD-oSC*SD,HIGH3SD=MD+oSC*SD,
                LOW95=MD-CI95,HIGH95=MD+CI95,
                TMBias_p=t.test(Y,mu=0)$p.val)

adboundr=cbind(adbound[,1:2],round(adbound[,3:13],2),TMBias_p=round(adbound[,14],4))
adboundNOr=cbind(adboundNO[,1:2],round(adboundNO[,3:13],2),TMBias_p=round(adboundNO[,14],4))
adboundr
adboundNOr

IAP=subset(AdErrorFile,Error_Category=="Imaging & Annotator Error"&Genus=="Porites sp.")
IAPab=subset(adbound,Error_Category=="Imaging & Annotator Error"&Genus=="Porites sp.")
plot(sort(IAP$TransitionMagnitude),ylim=c(-1000,1000))
abline(h=c(IAPab$LOW3SD,IAPab$HIGH3SD))

ggplot(AdErrorFile,aes(x=Y))+geom_histogram(binwidth=5)+
  geom_point(aes(y=1),data=subset(AdErrorFile,Outlier==T),color="red")+
  geom_point(aes(y=10),data=subset(AdErrorFile,AnnotatorNoteFlag==T),color="magenta")+
  facet_grid(Genus~Error_Category)


###Tables
library(sjPlot)
library(stargazer)
library(htmlTable)
visualoutput="C:/Users/Thomas.Oliver/WORK/Projects/Vital Rates/TechMemoVisuals/"
ErrSiteTab=ddply(AdErrorFile,.(Genus,Error_Category),summarize,N_Sites=length(unique(Site)),N_Colony=length(Y))
ErrSiteTab=ErrSiteTab[order(ErrSiteTab$Genus,ErrSiteTab$Error_Category),]

ft=flextable(data = ErrSiteTab)
ft = set_header_labels(ft,values=list(
  Error_Category="Error Category",
  Genus="Genus",
  N_Sites="N. Sites",
  N_Colony="N. Colonies"))
ft=merge_v(ft,j="Genus")
ft=autofit(ft)
#ft=theme_alafoli(ft)
#ft=border_inner_v(ft,border=fp_border(width=4),part = "header")
#ft=border_outer(ft,border=fp_border(width=4),part = "header")
ft

docx_file <- tempfile(pattern="ErrorSiteCount",tmpdir = visualoutput,fileext = ".docx")
save_as_docx("Error Count Table" = ft, path = docx_file)

adcol=c("Genus","Error_Category","N","MN","CI95","TMBias_p")

ftad=flextable(data = adboundNOr,col_keys = adcol)
ftad = set_header_labels(ftad,values=list(
  Error_Category="Error Category",
  Genus="Genus",
  N="N Colonies",
  MN=paste0("Mean ",YPlotLabel),
  CI95="+/- 95% CI",
  TMBias_p="Non-Zero Bias\n P-value"))
ftad=merge_v(ftad,j="Genus")
ftad=autofit(ftad)

docx_file <- tempfile(pattern=outname,
                      tmpdir = visualoutput,
                      fileext = ".docx")
save_as_docx(ftad, path = docx_file)



###
#Plot X,Y
pXY=ggplot(AdErrorFile)+
  geom_point(aes(x=log(StartingSize),
                 y=log(EndingSize),
                 color=Error_Category,
                 shape=Outlier))+
  geom_abline()+
  ylab("Size at T=1 (cm^2)")+
  xlab("Size at T=0 (cm^2)")+
  facet_grid(Error_Category~Genus)+
  scale_x_continuous(breaks=log(c(20,50,150,400,1000,3000)),labels=c(20,50,150,400,1000,3000))+
  scale_y_continuous(breaks=log(c(20,50,150,400,1000,3000)),labels=c(20,50,150,400,1000,3000))+
  scale_color_viridis(discrete = T)+
  scale_fill_viridis(discrete = T)+
  theme_bw()+scale_shape_manual(values=c(16,4))+
  ggtitle("Colony-Level Growth and Estimation Errors")+
  theme(legend.position = "none")
#pXY

#Plot Histograms
pHist=ggplot(AdErrorFileNO)+
  geom_histogram(aes(x=Y,
                 fill=Error_Category),bins=100)+
  geom_point(data=AEFYES,aes(x=Y),y=5,shape=4)+
  geom_vline(xintercept=0,color="black")+
  geom_rect(data=adboundNO,aes(xmin=LOW95,xmax=HIGH95,ymin=0,ymax=Inf),fill="lavender",color="gray25",alpha=0.75)+
  ylab("Count of Colonies")+
  xlab(YPlotLabel)+
  facet_grid(Error_Category~Genus,scales = "free_y")+
  theme_bw()+
  scale_fill_viridis(discrete = T)+
  scale_x_continuous(limits=c(min(AdErrorFileNO$Y),max(AdErrorFileNO$Y)))+
  ggtitle("Colony-Level Growth and Estimation Errors")+
  theme(legend.position = "none")
#pHist

#Plot Violin/Boxplots
pViBox=ggplot(AdErrorFileNO,aes(x=Error_Category,
                               y=Y,
                               fill=Error_Category))+
  geom_hline(yintercept = 0)+
  geom_violin(width=0.9,scale="width",alpha=0.33)+
  geom_boxplot(width=0.33,color="black",size=1)+
  geom_point(data=AEFYES,shape=4)+
  geom_text(aes(label=N,y=MD),data=adboundNO,color="white",size=4)+
  xlab("Error Category")+
  ylab(YPlotLabel)+
  facet_grid(.~Genus,scales = "free_y")+
  scale_fill_viridis(discrete = T)+
  theme_bw()+
  scale_y_continuous(limits=c(min(AdErrorFileNO$Y),max(AdErrorFileNO$Y)))+
  ggtitle("Colony-Level Growth and Estimation Errors")+
  theme(legend.position = "none")
#pViBox


#Plot Just Boxplots
pBox=ggplot(AdErrorFileNO,aes(x=Error_Category,
                                y=Y,
                                fill=Error_Category))+
  geom_hline(yintercept = 0)+
  geom_boxplot(width=0.33,color="black",size=1)+
  geom_point(data=AEFYES,shape=4)+
  geom_text(aes(label=N,y=MD),data=adboundNO,color="white",size=4)+
  xlab("Error Category")+
  ylab(YPlotLabel)+
  facet_grid(.~Genus,scales = "free_y")+
  scale_fill_viridis(discrete = T)+
  theme_bw()+
  scale_y_continuous(limits=c(min(AdErrorFileNO$Y),max(AdErrorFileNO$Y)))+
  ggtitle("Colony-Level Growth and Estimation Errors")+
  theme(legend.position = "none")
#pBox

sc=.85
ggsave(plot = pXY,
       filename = paste0("C:/Users/Thomas.Oliver/WORK/Projects/Vital Rates/ErrorRates/ErrorXY_",outname,".jpg"),
       width=sc*8.5,height=sc*11)
ggsave(plot = pHist,
       filename = paste0("C:/Users/Thomas.Oliver/WORK/Projects/Vital Rates/ErrorRates/ErrorHist_",outname,".jpg"),
       width=sc*8.5,height=sc*11)
sc=1.25
ggsave(plot = pViBox,
       filename = paste0("C:/Users/Thomas.Oliver/WORK/Projects/Vital Rates/ErrorRates/ErrorViBox_",outname,".jpg"),
       width=sc*11,height=sc*8.5)
ggsave(plot = pBox,
       filename = paste0("C:/Users/Thomas.Oliver/WORK/Projects/Vital Rates/ErrorRates/ErrorBox_",outname,".jpg"),
       width=sc*11,height=sc*8.5)


#Pct Change with Size
SizevsPct=ggplot(ErrorFile,aes(x=StartingSize,y=TransitionMagnitude,color=Error_Category,fill=Error_Category))+
  geom_point(size=1,data=AdErrorFile)+
  geom_point(size=.25,color="gray")+
  stat_smooth(method="loess",span=1.5)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = c(19.6,50.3),color="black",lty=2)+
  scale_y_continuous(limits=c(-100,150))+
  scale_x_log10(breaks=c(3,20,50,100,200,500,1000),limits=c(3,1000))+
  scale_color_viridis(discrete = T)+
  scale_fill_viridis(discrete = T)+
  ggtitle("Percentage Size Change as Function of Coral Patch Size")+
  xlab("Colony Size at T=0 (cm^2)")+
  ylab("Colony Percentage Size Change (% cm^2)")+
  annotate(geom="text",x=7,y=-75,label="'Juveniles'")+
  annotate(geom="text",x=100,y=-75,label="'Adults'")+
  facet_grid(Error_Category~Genus)+theme_bw()+theme(legend.position = "none")
sc=0.85
ggsave(plot = SizevsPct,
       filename = paste0("C:/Users/Thomas.Oliver/WORK/Projects/Vital Rates/ErrorRates/ErrorSizeByPct.jpg"),
       width=sc*8.5,height=sc*11)

# Other Plots -------------------------------------------------------------



ggplot(AdErrorFileNO)+
  geom_point(aes(x=log(StartingSize),
                 y=Y,
                 color=Error_Category,
                 shape=Outlier))+
  geom_hline(yintercept=0)+
  geom_hline(data=lines,aes(yintercept = MN+CI95),color="red",lty=2)+
  geom_hline(data=lines,aes(yintercept = MN),color="red",lty=1)+
  geom_hline(data=lines,aes(yintercept = MN-CI95),color="red",lty=2)+
  ylim(c(-200,200))+
  ylab(YPlotLabel)+
  xlab("Size at T=0 (log-cm^2)")+
  facet_grid(Error_Category~Genus,scales="free")+
  theme_bw()+scale_shape_manual(values=c(16,4))
#Is there bias? 


ggplot(AdErrorFile)+
  geom_histogram(aes(x=log10(1+StartingSize)),bins=100,col="red")+
  geom_histogram(aes(x=log10(1+EndingSize)),bins=100,col="blue")+
  geom_vline(xintercept=0)+
  facet_wrap(Error_Category~Genus_Code)+theme_bw()

ggplot(AdErrorOnly,
       aes(StartingSize,(Log2Ratio),color=Outlier))+
  geom_point()+
  facet_grid(Error_Category~Genus)+theme_bw()

# steps=c(1,5,10,20,50,seq(100,3000,by=100))
# DF=NULL
# for(i in 1:length(steps)){
#   DF=rbind(DF,table(cut(ErrorFile$Log2Ratio,breaks=c(-Inf,-steps[i],steps[i],Inf))))
# }
# DF=cbind(steps,DF)

plot(DF[,1],DF[,2],type="b",ylim=c(-10,10))
points(DF[,1],-DF[,4],type="b")


ErrorFile$Y=abs(ErrorFile$Log2Ratio)
lines=ddply(ErrorFile,.(Genus,Error_Category),summarize,
            MN=median(Y),CI95=1.96*sd(Y)/sqrt(length(Y)),LOW=MN-CI95,HIGH=MN+CI95)
ggplot(AdErrorFile,aes(y=(Log2Ratio)))+
  geom_boxplot(width=.25)+
  geom_hline(data=lines,aes(yintercept = MN+CI95),color="red")+
  geom_hline(data=lines,aes(yintercept = MN-CI95),color="red")+
  ylim(c(-2,2))+
  geom_hline(yintercept =0)+
  facet_grid(Genus~Error_Category)


ErrorFile$Y=ErrorFile$DiffOfLog
lines=ddply(ErrorFile,.(DataOrError,Genus_Code),summarize,MN=median(Y),CI95=1.96*sd(Y)/sqrt(length(Y)),LOW=MN-CI95,HIGH=MN+CI95)
ggplot(subset(ErrorFile,Genus_Code!="MOSP"))+
  geom_histogram(aes(x=DiffOfLog,fill=DataOrError),bins=100,color="black")+
  geom_vline(xintercept = 0,color="black")+
  geom_vline(data=subset(lines,Genus_Code!="MOSP"),aes(xintercept = LOW),color="red")+
  geom_vline(data=subset(lines,Genus_Code!="MOSP"),aes(xintercept = HIGH),color="red")+
  facet_grid(DataOrError~Genus_Code,scales="free_y")+theme_bw()


ggplot(ErrorFile,aes(StartingSize,EndingSize))+
  geom_point()+
  geom_abline()+
  facet_grid(DataOrError~Genus_Code)


ggplot(ErrorFile,aes(x=DataOrError,y=TransitionMagnitude))+
  geom_boxplot()+
  ylim(c(-20,20))+
  geom_hline(yintercept =0)+
  facet_grid("Genus_Code")


ggplot(ErrorFile,aes(x=DataOrError,y=Log2Ratio))+
  geom_boxplot()+
  ylim(c(-2,2))+
  geom_hline(yintercept =0)+
  facet_grid("Genus_Code")


ggplot(ErrPlot,aes(x=DataOrError,y=Log2Ratio_Change))+geom_violin(binwidth=5)+facet_wrap("Genus_Code")
bounds=ddply(VGS_GSerr,.(Genus_Code),summarize,low95=quantile(PercentChange,.05),high95=quantile(PercentChange,.95))

ggplot(VGR_GSerr,aes(x=PercentChange))+
  geom_histogram()+
  facet_grid(Genus_Code~Annotator_Head)

VRGde_nffmr=subset(VitalRate_Growth,TransitionType%in%c("GROWTH","SHRINK"))

VRGd=subset(VitalRate_Growth,DataOrError=="DATA")
MaxDiamBreaks=c(0,5,10,20,30,50,100)
VRGd$SizeBin=cut(VRGd$MaxDiam,breaks=MaxDiamBreaks,include.lowest=T)#,labels=c("A.Small","B.M1","C.M2","D.Large"))
VRGd$Recruit=as.numeric(VRGd$TransitionType=="RECR")
VRGd$Mortality=as.numeric(VRGd$TransitionType=="MORT")
VRGd_nr=subset(VRGd,TransitionType%in%c("GROWTH","SHRINK","FISSION","FUSION","FUSION_FISSION","MORT"))
VRGd_nmr=subset(VRGd,TransitionType%in%c("GROWTH","SHRINK","FISSION","FUSION","FUSION_FISSION"))
VRGd_nff=subset(VRGd,TransitionType%in%c("GROWTH","SHRINK","MORT","RECR"))
VRGd_nffr=subset(VRGd,TransitionType%in%c("GROWTH","SHRINK","MORT"))
VRGd_nffrm=subset(VRGd,TransitionType%in%c("GROWTH","SHRINK"))

# #Plot Colony Graphs ---------------------------------------------------------------
scale=250
jpeg(filename = paste0(basepath,"Figures/OAH010_Allgraphs.jpg"),width = 4*scale,height=3*scale)
par(mfrow=c(1,2))
plot.colgraph(subg_Eind(gOAH010_POCS,
                        edge_indices = sample(1:length(E(gOAH010_POCS)),16)),
              layoutfun = layout_nicely)
plot.colgraph(gOAH010_POSP)
dev.off()

scale=250
jpeg(filename = paste0(basepath,"Figures/D2_PMEA_Error_All36graphs.jpg"),width = 4*scale,height=3*scale)
par(mfrow=c(1,2))
plot.colgraph(gD2_POCS)
plot.colgraph(gD2_POSP)
dev.off()


# #Plot SUB Graphs ---------------------------------------------------------------
subv=sort(c("P_054","P_053","P_006","P_033","P_035","P_003","P_002","P_001"))
# nv=length(V(gOAH010))
subgOAH010=subg_Vname(gOAH010,subv)
scale=250
jpeg(filename = paste0(basepath,"Figures/OAH010_PMEA_Data_SUBgraphs.jpg"),width = 4*scale,height=3*scale)
plot.igraph(subgOAH010,layout=layout_as_tree(subgOAH010),
            edge.width=SizeFunc(E(subgOAH010)$AreaChange_cm2,lowest=.05,highest=3),
            edge.arrow.width=SizeFunc(abs(E(subgOAH010)$AreaChange_cm2),lowest=.05,highest=3),
            edge.color=EdgeColor[E(subgOAH010)$TRANSITION_TYPE],
            vertex.color=NodeColor[V(subgOAH010)$Type],
            vertex.size=SizeFunc(V(subgOAH010)$Area),
            vertex.shape=c("circle","square")[V(subgOAH010)$TimePt],
            vertex.label.cex=0.75,
            vertex.label.color="black",
            vertex.label.dist=1.5)
dev.off()

# Plot Growth metrics vs Error ####
#Growth vs Error - OAH010, D2
ErrorBounds_POCS=quantile(subset(VRGde_nffmr,DataOrError=="ERROR"&GenusCode=="POCS")$TransitionMagnitude,c(.05,.95))
ErrorBounds_POSP=quantile(subset(VRGde_nffmr,DataOrError=="ERROR"&GenusCode=="POSP")$TransitionMagnitude,c(.05,.95))
DvEraw_POCS=ggplot(subset(VRGde_nffmr,GenusCode=="POCS"),
                   aes(y=TransitionMagnitude,x=DataOrError,color=DataOrError))+
  geom_boxplot()+
  geom_jitter(height = 0,width=.1)+
  geom_hline(yintercept = ErrorBounds_POCS,color="red")+
  ylab("Measured Change in Planar Area (cm2)")+xlab("")+
  ggtitle(paste0("POCS: Raw Measured Growth vs Measured Error\n No Fusion,Fission,Mortality,Recruitment\nError Bounds: ",
                 paste0(round(ErrorBounds_POCS,2),collapse=" - ")," cm^2"))+
  #scale_color_aaas()+
  ylim(c(-300,300))+theme_bw()
DvEraw_POSP=ggplot(subset(VRGde_nffmr,GenusCode=="POSP"),
                   aes(y=TransitionMagnitude,x=DataOrError,color=DataOrError))+
  geom_boxplot()+
  geom_jitter(height = 0,width=.1)+
  geom_hline(yintercept = ErrorBounds_POSP,color="red")+
  ylab("Measured Change in Planar Area (cm2)")+xlab("")+
  ggtitle(paste0("POSP: Raw Measured Growth vs Measured Error\n No Fusion,Fission,Mortality,Recruitment\nError Bounds: ",
                 paste0(round(ErrorBounds_POSP,2),collapse=" - ")," cm^2"))+
  #scale_color_aaas()+
  ylim(c(-300,300))+theme_bw()
DvEraw=ggarrange(DvEraw_POCS,DvEraw_POSP)
scale=.7
ggsave(filename = paste0(basepath,"Figures/Growth_DataVsError_v2.png"),
       plot=DvEraw,width=16*scale,height=scale*9)

# Plot Transition Type By Size ####
#Transition Type By Size
TTbySc=ggplot(VRGd,aes(x=MaxDiam,fill=TransitionType))+
  geom_histogram(bins=20)+
  facet_grid("GenusCode")+
  scale_fill_brewer(type = "qual")+theme_bw()
ggsave(filename = paste0(basepath,"Figures/TransitionTypeBySize_Cont_v2.png"),plot=TTbySc,width=4*scale,height=scale*3)


# Plot Growth T, vs T+1 ####
#GRowth T, vs T+1
gT0_T1=ggplot(subset(VRGd,Genus_Code%in%c("POCS","POSP")),
              aes(x=StartingSize,y=EndingSize,color=TransitionType,shape=TransitionType,group=ColonyID))+
  geom_point()+
  geom_line(alpha=.15,color="black")+
  geom_abline(slope=1)+
  scale_x_sqrt(breaks=c(10,50,100,250,500,1000,2000))+
  scale_y_sqrt(breaks=c(10,50,100,250,500,1000,2000))+
  scale_shape_manual(values=c(3,4,1,0,6,2,11))+
  facet_wrap(c("Genus_Code","Site"),scales="free",nrow = 2)+
  theme_bw()+
  xlab("Colony Size at Time T (cm^2)")+
  ylab("Colony Size at Time T+1 (cm^2)")

scale=2
ggsave(filename = paste0(basepath,"Figures/TransitionTypeSizeTSizeTp1_v3.png"),
       plot=gT0_T1,width=6*scale,height=scale*3)

gT0_T1z=ggplot(subset(VRGd,Genus_Code%in%c("POCS","POSP","MOSP")),
               aes(x=StartingSize,y=EndingSize,color=TransitionType,shape=TransitionType))+
  geom_point()+
  geom_abline(slope=1)+
  scale_x_sqrt(breaks=c(10,50,100,250,500,1000,2000),limits=c(0,400))+
  scale_y_sqrt(breaks=c(10,50,100,250,500,1000,2000),limits=c(0,400))+
  scale_shape_manual(values=c(3,4,1,0,6,2,11))+
  facet_grid(c("Genus_Code","Site"))+theme_bw()+
  xlab("Colony Size at Time T (cm^2)")+
  ylab("Colony Size at Time T+1 (cm^2)")
scale=3
ggsave(filename = paste0(basepath,"Figures/TransitionTypeSizeTSizeTp1_ZOOMv3.png"),
       plot=gT0_T1z,width=6*scale,height=scale*3)

K01 <- subset(VRGd_nff, Site == "MAI_SIO_K01")
k01plot <- ggplot(subset(K01,Genus_Code%in%c("POSP","MOSP")),
                          aes(x=StartingSize,y=EndingSize,color=TransitionType,shape=TransitionType))+
  geom_point()+
  geom_abline(slope=1)+
  scale_x_sqrt(breaks=c(10,50,100,250,500,1000,2000),limits=c(0,400))+
  scale_y_sqrt(breaks=c(10,50,100,250,500,1000,2000),limits=c(0,400))+
  scale_shape_manual(values=c(3,4,1,0,6,2,11))+
  scale_color_manual(values = c("#009E73","#FF0033","#56B4E9","#FF9900"))+
  facet_grid(c("Genus_Code","StartingDate"))+theme_bw()+
  xlab("Colony Size at Time T (cm^2)")+
  ylab("Colony Size at Time T+1 (cm^2)")


# Plot Recruit Rate ####
RecAdDen=ddply(VRGd,.(Site,Genus_Code),summarize,
         Nrecruit=sum(Recruit),
         Nrec_area=Nrecruit/mean(AreaSurveyed,na.rm=T),
         Nad=length(unique(ColonyID))-Nrecruit,
         Nad_area=Nad/mean(AreaSurveyed,na.rm=T),
         Ntp=1,
         Nrec_Nad=(Nrec_area/Nad_area)/Ntp)

RecDen=ggplot(RecAdDen,aes(x=Site,y=Nrec_area,fill=Genus_Code))+
  geom_hline(yintercept = 0,color="black")+
  geom_bar(stat="Identity",position="dodge",width=.5)+
  xlab("Site")+
  ylab("Recruits per Area (N/m^2)")+
  ggtitle("Site/Genus Specific Recruit Density")+theme_bw()
ggsave(filename = paste0(basepath,"Figures/RecruitDensity_v3.png"),plot=RecDen,width=4*scale,height=scale*3)

RecOverAdDen=ggplot(RecAdDen,aes(x=Site,y=Nrec_Nad,fill=Genus_Code))+
  geom_hline(yintercept = 0,color="black")+
  geom_bar(stat="Identity",position="dodge",width=.5)+
  xlab("Site")+
  ylab("Recruits/Adult Colony Ratio")+
  ggtitle("Site/Genus Specific Recruits per Adult Colony")
ggsave(filename = paste0(basepath,"Figures/RecruitAdultRatio.png"),plot=RecOverAdDen,width=4*scale,height=scale*3)

# Plot Size Specific Mortality Rate ####
modMort=glm(formula=Mortality~StartingSize:Genus_Code:Site,data=VRGd_nffr,family=binomial)
summary(modMort)
#Prediction for each Site, Genus, Size Combo
SizeRange=seq(0,1500,length.out=100)
uGenera=unique(VRGd_nffr$Genus_Code)
uSites=unique(VRGd_nffr$Site)
nSizeReps=length(uGenera)*length(uSites)
nGenReps=length(SizeRange)*length(uSites)
nSiteReps=length(uGenera)*length(SizeRange)
pred_VRGd_nffr=data.frame(Site=sort(rep(uSites,nSiteReps)),
                          Genus_Code=rep(sort(rep(uGenera,length(SizeRange))),length(uSites)),
                          StartingSize=rep(SizeRange,nGenReps))
#subset predictions within size range of data...
uS=unique(pred_VRGd_nffr$Site)
uG=unique(pred_VRGd_nffr$Genus_Code)
predsub=NULL
SR=ddply(VRGd_nffr,.(Site,Genus_Code),summarize,
         Size_min=quantile(StartingSize,0,na.rm=T),
         Size_q01=quantile(StartingSize,.01,na.rm=T),
         Size_q99=quantile(StartingSize,.99,na.rm=T),
         Size_max=quantile(StartingSize,1,na.rm=T))

for(i in 1:length(uS)){
  for(j in 1:length(uG)){
    this_SG=subset(pred_VRGd_nffr,
                   Genus_Code==uG[j]&
                     Site==uS[i]&
                     StartingSize>subset(SR,Genus_Code==uG[j]&Site==uS[i])$Size_q01&
                     StartingSize<subset(SR,Genus_Code==uG[j]&Site==uS[i])$Size_q99)
    predsub=rbind(predsub,this_SG)
  }
}


pModMort=predict(modMort,se=T,
                 newdata=predsub,
                 type="response")
predsub$Mortality=pModMort$fit
predsub$Mort_min=pModMort$fit-pModMort$se.fit
predsub$Mort_max=pModMort$fit+pModMort$se.fit
SSM=ggplot()+
  geom_vline(xintercept=pi*(2.5^2),color="gray50")+
  geom_hline(yintercept=c(0,1),color="gray50")+
  geom_jitter(data=VRGd_nffr,aes(x=StartingSize,y=Mortality,color=Site),height=.05)+
  geom_ribbon(data=predsub,aes(x=StartingSize,ymin=Mort_min,ymax=Mort_max,fill=Site),alpha=.1)+
  geom_line(data=predsub,aes(x=StartingSize,y=Mortality,color=Site),size=1.25,lty=1)+
  #  geom_line(data=predsub,aes(x=StartingSize,y=Mortality,color=Site),size=1,lty=1)+
  ylab("Mortality (Prob.)")+
  xlab("Starting Colony Area (cm2)")+
  ggtitle("Size-Dependent Mortality")+
  facet_grid(.~GenusCode)+
  scale_x_sqrt(breaks=c(0,1,25,100,200,500,1000),limits=c(0,1000))+theme_bw()
scale=1.5
ggsave(filename = paste0(basepath,"Figures/Size-Specific Mortality_v2.png"),
       plot=SSM,width=6*scale,height=scale*3)

#Size Specific Growth Rate
modGrow=glm(formula=TransitionRate~StartingSize:GenusCode:Site,data=VRGd_nffrm,family=gaussian)
summary(modGrow)
modSSG=glm(formula=SizeSpecificRate~StartingSize:GenusCode:Site,data=VRGd_nffrm,family=gaussian)
summary(modSSG)
SizeRange=seq(0,1000,length.out=100)
pModGrow=predict(modGrow,se=T,
                 newdata=predsub,
                 type="response")
pModSSG=predict(modSSG,se=T,
                newdata=predsub,
                type="response")
predsub$Growth=pModGrow$fit
predsub$Growth_min=pModGrow$fit-pModGrow$se.fit
predsub$Growth_max=pModGrow$fit+pModGrow$se.fit
predsub$SSG=pModSSG$fit
predsub$SSG_min=pModSSG$fit-pModSSG$se.fit
predsub$SSG_max=pModSSG$fit+pModSSG$se.fit


GbyS=ggplot()+
  geom_jitter(data=VRGd_nffrm,aes(x=StartingSize,y=TransitionRate,color=Site),size=1,width=1,height=1)+
  geom_ribbon(data=predsub,aes(x=StartingSize,ymin=Growth_min,ymax=Growth_max,fill=Site),alpha=.1)+
  geom_line(data=predsub,aes(x=StartingSize,y=Growth,color=Site))+
  ylab("Growth Rate (cm2/yr)")+
  xlab("Colony Area (cm2)")+
  ggtitle("Size Dependent Growth")+
  geom_hline(yintercept=0,color="black")+
  geom_vline(xintercept=pi*(2.5^2),color="gray50")+
  facet_grid(.~GenusCode)+
  scale_y_continuous(limits=c(-60,100))+
  scale_x_sqrt(breaks=c(1,25,100,200,500,1000))+theme_bw()
scale=1.5
ggsave(filename = paste0(basepath,"Figures/Growth By Size_v2.png"),
       plot=GbyS,width=6*scale,height=scale*3)


SSG=ggplot()+
  geom_jitter(data=VRGd_nffrm,aes(x=StartingSize,y=SizeSpecificRate,color=Site),size=1,width=1,height=1)+
  geom_ribbon(data=predsub,aes(x=StartingSize,ymin=SSG_min,ymax=SSG_max,fill=Site),alpha=.1)+
  geom_line(data=predsub,aes(x=StartingSize,y=SSG,color=Site))+
  ylab("Realtive Growth Rate (% Growth)")+
  xlab("Colony Area (cm2)")+
  ggtitle("Relative Growth by Size")+
  geom_hline(yintercept=0,color="black")+
  geom_vline(xintercept=pi*(2.5^2),color="gray50")+
  facet_grid("GenusCode")+
  #  scale_y_continuous(limits=c(-300,300))+
  scale_x_sqrt(breaks=c(1,25,100,200,500,1000))+theme_bw()
ggsave(filename = paste0(basepath,"Figures/SizeSpecificGrowth_v2.png"),plot=SSG,width=4*scale,height=scale*3)


##
vrO10=subset(VRGd_nffrm,Site=="OAH010"&GenusCode=="POSP")
vrO10$MidpointSize=((vrO10$StartingSize+vrO10$EndingSize)/2)
vrO10$SizeSpecificRate_mid=vrO10$TransitionRate/vrO10$MidpointSize

plot(sqrt(vrO10$MidpointSize),vrO10$TransitionRate)

table(VRGd_nffrm$GenusCode,VRGd_nffrm$Site)

# IPMpack PlayBall --------------------------------------------------------
library(IPMpack)
head(VRGd)
Select_Genus="POSP"
Select_Sites=c("OAH010")
PlotTitle=paste0(Select_Genus," at ",paste(Select_Sites,collapse = ","))
TargetSite_Taxon=subset(VRGd,GenusCode==Select_Genus&Site%in%Select_Sites)
MinFecSize=80
FEC_METHOD="FLAT"#"FLAT" or "SAMPLED"

# BUild Dataframe ####
SiteTaxTrans=TargetSite_Taxon[,c("Site","GenusCode","T0_PatchName","T1_PatchName","StartingSize","EndingSize","Mortality","TransitionType")]
SiteTaxTrans$Site_Gen=paste0(SiteTaxTrans$Site,"_",SiteTaxTrans$GenusCode)
#Basics, T,T+1 Size and Survival (1/0)
STipm=SiteTaxTrans[,c("StartingSize","EndingSize","Mortality")]
names(STipm)=c("size","sizeNext","surv")
STipm$surv=as.numeric(!STipm$surv)


#TaxonSpecific RecruitsPerSize
#basically, 0 until 20 cm squared, and then a falt value of site-level Recruit/Adult colony ratio
Site_Gen_RArat=RecAdDen$Nrec_Nad
names(Site_Gen_RArat)=paste0(RecAdDen$Site,"_",RecAdDen$GenusCode)
STipm$fec1=0
if(FEC_METHOD=="SAMPLED"){
  fec2_mean=Site_Gen_RArat[SiteTaxTrans$Site_Gen]
  fec2_sd=Site_Gen_RArat[SiteTaxTrans$Site_Gen]*.10
  fec2_draw=rep(NA,length(fec2_mean))
  for(i in 1:length(fec2_draw)){fec2_draw[i]=rnorm(n = 1,mean = fec2_mean[i],sd = fec2_sd[i])}
  STipm$fec2=fec2_draw 
}else{
  STipm$fec2=Site_Gen_RArat[SiteTaxTrans$Site_Gen]
}
STipm$fec1[STipm$size>=MinFecSize]=1
STipm$fec2[STipm$fec1==0]=0

#Stage/StageNext/OffspringNext For Rec/Mort/Fission/Fusion Cases
STipm$stage="continuous"
STipm$stageNext="continuous"
STipm$offspringNext=NA

#Now begin blocking off Rec/Mort
#Rec
STipm$size[which(SiteTaxTrans$TransitionType=="RECR")]=NA
#STipm$stage[which(SiteTaxTrans$TransitionType=="RECR")]="larvae"
STipm$offspringNext[which(SiteTaxTrans$TransitionType=="RECR")]="sexual"
#Mort
STipm$sizeNext[which(SiteTaxTrans$TransitionType=="MORT")]=NA
STipm$stageNext[which(SiteTaxTrans$TransitionType=="MORT")]="dead"

#Now Fission/Fusion
#Fission, group by head node (ie starting patch)
uT0P=unique(SiteTaxTrans$T0_PatchName[SiteTaxTrans$TransitionType%in%c("FISSION","FUSION_FISSON")])
#for each head node with fission as an transition...
for(i in 1:length(uT0P)){
  #Find all in STipm space
  fis_i=which(SiteTaxTrans$T0_PatchName==uT0P[i])
  #Find biggest piece in STipm space
  fis_i_main=fis_i[which.max(SiteTaxTrans$EndingSize[fis_i])]
  fis_i_minor=setdiff(fis_i,fis_i_main)
  #update size, and offspringNext
  #STipm$size[fis_i_minor]=NA
  STipm$offspringNext[fis_i_minor]="clonal"
}
#Fusion, group by head node (ie starting patch)
uT1P=unique(SiteTaxTrans$T1_PatchName[SiteTaxTrans$TransitionType%in%c("FUSION","FUSION_FISSON")])
#for each head node with fission as an transition...
for(i in 1:length(uT1P)){
  #Find all in STipm space
  fus_i=which(SiteTaxTrans$T1_PatchName==uT1P[i])
  #Find biggest piece in STipm space
  fus_i_main=fus_i[which.max(SiteTaxTrans$StartingSize[fus_i])]
  fus_i_minor=setdiff(fus_i,fus_i_main)
  #update size, and offspringNext
  STipm$sizeNext[fus_i_minor]=NA
  STipm$surv[fus_i_minor]=0
  STipm$stageNext[fus_i_minor]="Fusion_Mortality"
}
#Reapply fec to new sizes
STipm$fec1[is.na(STipm$size)]=0
STipm$fec2[is.na(STipm$size)]=0
STipm$stage=factor(STipm$stage)
STipm$stageNext=factor(STipm$stageNext)
STipm$offspringNext=factor(STipm$offspringNext)


#Given stable Data, Run Model Objects ####

#Growth/Surv
gr1=makeGrowthObj(dataf=STipm,Formula=sizeNext~size+size2)
sv1=makeSurvObj(dataf = STipm,Formula=surv~size+size2)

Pmatrix=makeIPMPmatrix(nBigMatrix = 50,
                       minSize=-5,maxSize=ceiling(1.1*max(STipm$sizeNext,na.rm=T)),
                       growObj = gr1,survObj = sv1,
                       correction="constant")
LE=meanLifeExpect(Pmatrix)
pTime=passageTime(mean(STipm$size,na.rm=T),Pmatrix)


#Sexual Vs Clonal Repro
table(STipm$offspringNext)
fs1=makeFecObj(dataf = STipm,Formula=fec2~size+size2)
fc1=makeClonalObj(dataf = STipm,Formula=sizeNext~size+size2)

Fmatrix=makeIPMFmatrix(nBigMatrix=50,minSize=-5,
                       maxSize = ceiling(1.1*max(STipm$sizeNext,na.rm=T)),
                       fecObj = fs1,correction = "constant")
Cmatrix=makeIPMCmatrix(nBigMatrix=50,minSize=-5,
                       maxSize = ceiling(1.1*max(STipm$sizeNext,na.rm=T)),
                       clonalObj = fc1,correction = "constant")

IPM=Pmatrix+Fmatrix#+Cmatrix
library(fields)


par(mfrow=c(1,2),bty="l",pty="m")
p1=picGrow(STipm,gr1)
p2=picSurv(STipm,sv1,ncuts=30)

par(mfrow=c(2,2),mar=c(4,4,6,2))
image.plot(t(Pmatrix),col=topo.colors(20),
           main=paste0(PlotTitle,"\nPmatrix: Growth and Survival"),xlab="T",ylab="T+1")
abline(a=0,b=1,col="gray50",lwd=2)
image.plot(t(Fmatrix),col=topo.colors(20),
           main=paste0(PlotTitle,"\nFmatrix: Sexual Fecundity"),xlab="T",ylab="T+1")
abline(a=0,b=1,col="gray50",lwd=2)
image.plot(t(Cmatrix),col=topo.colors(20),
           main=paste0(PlotTitle,"\nCmatrix: Clonal Fecundity"),xlab="T",ylab="T+1")
abline(a=0,b=1,col="gray50",lwd=2)
image.plot(t(IPM),col=topo.colors(20),
           main=paste0(PlotTitle,"\nIPM Matrix (No Clonal): Lambda=",round(Re(eigen(IPM)$value[1]),4)),xlab="T",ylab="T+1")
abline(a=0,b=1,col="gray50",lwd=2)

sensitivity=sens(IPM)
elasticity=elas(IPM)

res=sensParams(growObj = gr1,survObj = sv1,fecObj = fs1,
               nBigMatrix = 50,minSize =-5,maxSize = max(STipm$sizeNext,na.rm=T))
res



diagnosticsPmatrix(Pmatrix,growObj = gr1,survObj = sv1)
plot(Pmatrix@meshpoints,LE)
plot(Pmatrix@meshpoints,pTime)

par(mfrow=c(1,1),bty="l",pty="m",mar=c(10,4,3,1))
#barplot(res$sens,main=expression("Parameter Sensitivity of "*lambda),las=2,cex.names=0.5)
barplot(res$elas,main=expression("Parameter Elasticity of "*lambda),las=2,cex.names=0.65)
