#IPM Analyses

rm(list=ls())
# Loading Libraries -------------------------------------------------------
library(tidyverse)
library(lubridate)
library(stringr) 
library(pROC)
library(ggpubr)
library(sf)
library(rgeos)
library(sp)
library(patchwork)
library(gridExtra)


# Functions -------------------------------------------------------
inv.logit <- function(x) {
  exp(x)/(1+exp(x))
}

################################################################################
###### The growth model ########################################################
################################################################################

growth_function <- function(y, x, g.int, g.slp, g.var) {
  mu <- (g.int + g.slp * x)
  sig <- sqrt(g.var)
  dnorm(y, mean=mu, sd=sig)
}

################################################################################
###### The survival model ######################################################
################################################################################


survival_function <- function(x, s.int, s.slp,interval_years) {
  #, s.int, s.slp)
  Line <- s.slp * x + s.int
  u <- inv.logit(Line)#^(1/interval_years) # From binomial glm, link = logit function
  return(u)
}

################################################################################
###### The reproduction model ##################################################
################################################################################

reproduction_function <- function(y, x, rec, rec.size=rec.size,ind = FALSE) {
  # 2*sqrt((10^-2.5)/pi)  # approx 6 cm  diameter coral
  if (ind == TRUE){
    out <- rep(0,length(x)) #changed
    out[1:4] <- rec/4} #changed
  
  #out[x < rec.size | y >= rec.size] <- 0}
  
  if (ind == FALSE){
    out <-  (10^x) * rec
    out[x < rec.size | y >= rec.size] <- 0 #if x is below recruitment size then it doesn't count 
  }
  return(out)
}


reproduction_function_josh <- function(y, x,rec, rec.size=rec.size) {
  out <-  (10^x) * rec
  out[x < rec.size | y >= rec.size] <- 0
  return(out)
}

################################################################################
########## ... Running the IPM ... #############################################
################################################################################

bigmatrix <- function(ModelParameters,MatrixVal) {
  rec <- ModelParameters$rec
  rec.size <- MatrixVal$rec.size
  g.int <- ModelParameters$g.int
  g.slp <- ModelParameters$g.slp
  g.var <- ModelParameters$g.var
  s.int <- ModelParameters$s.int
  s.slp <- ModelParameters$s.slp
  #iy<- ModelParameters$Interval_Years
  
  #Reproduction_kernel <- delta_size * outer(y, y, 
  #reproduction_function, rec)
  Reproduction_kernel <- MatrixVal$delta_size * outer(MatrixVal$y, MatrixVal$y, 
                                            reproduction_function,rec=rec, 
                                            rec.size=MatrixVal$rec.size)
  
  Growth_kernel <- MatrixVal$delta_size * outer(MatrixVal$y, MatrixVal$y, growth_function, g.int=g.int, 
                                      g.slp=g.slp, g.var=g.var)
  ## Growth Kernel
  # only works inside the function or if you have done this by had
  #image(y,y,t(Growth_kernel))
  
  Survival_kernel <- survival_function(MatrixVal$y, s.int=s.int, s.slp=s.slp)#,interval_years = iy)
  
  P <- Growth_kernel
  i <- 1:MatrixVal$n
  
  P[,i] <- Growth_kernel[,i]*Survival_kernel[i]
  
  K <- P + Reproduction_kernel
  
  lam <- Re(eigen(K)$values[1])
  w <- Re(eigen(K)$vectors[,1])
  stabledist <- w/sum(w) # right eigenvector/stable distribution
  v.eigen <- Re(eigen(t(K))$vectors[,1]) #left eigenvector
  reprodvalue <- v.eigen/v.eigen[1] #reproductive value
  
  #compute sensitivity and elasticity matrices using the eigenvectors and eigenvalue
  v.dot.w <- sum(stabledist * reprodvalue) * MatrixVal$delta_size
  sens <- outer(reprodvalue,stabledist) / v.dot.w
  elas <- matrix(as.vector(sens) * as.vector(K) / lam, nrow = MatrixVal$n)
  
  return(list(K=K, Gk=Growth_kernel,Sk=Survival_kernel,Pk=P,lam=lam, w=stabledist, v=reprodvalue,
              v.dot.w=v.dot.w, sens=sens, elas=elas))
}

# Loading Data -------------------------------------------------------
Name="HA_MA_AS_Models"
loadnames3=load("./Data/ColonyTransitions/Script_Step3_DataPackage.rdata");loadnames3
loadnameSIG=load( file = sprintf("./Data/ModelData/%s_allmodfits_noNAs.rdata",Name));loadnameSIG
loadnameREG=load( file = sprintf("./Data/ModelData/%s_RegionalModFits_noBADs.rdata",Name));loadnameREG

MP=full_join(ModelParams_Regional,ModelParams_SIG)


# Run Models  -------------------------------------------------------

#First Pull Together Models from Current and Former Work
# Global mesh variables 
MatrixVal=list(
  #min.size <- 0.9*min(c(ColTrans$log10_SS, ColTrans$log10_ES), na.rm = T)
  min.size = -0.5,   # 2*sqrt((10^-2.5)/pi)  # approx 6 cm  diameter coral
  #max.size = 2,     # 2*sqrt((10^2)/pi)     # approx 11 m diameter coral (ridiculous, but avoids boundary issues)
  max.size = 1.1*max(c(ColTrans$log10_SS, ColTrans$log10_ES), na.rm = T),
  # max ~ 6.4 m colony diameter
  n = 50, #mesh size / number of cells in the discretized kernel
  rec.size = -0.1)#,  # 2*sqrt((10^-0.1)/pi)  # approx 1 cm  diameter coral #2*sqrt((10^1.29303)/pi)  # approx 5 cm  diameter coral
  
MatrixVal$bin_size = MatrixVal$min.size + c(0:MatrixVal$n) * (MatrixVal$max.size - MatrixVal$min.size)/MatrixVal$n #boundary points (the edges of the cells defining the kernel 
MatrixVal$y = 0.5 * (MatrixVal$bin_size[1:MatrixVal$n]+MatrixVal$bin_size[2:(MatrixVal$n+1)]) #mesh points (midpoints of cells)
MatrixVal$I = MatrixVal$y >= MatrixVal$rec.size
MatrixVal$delta_size = MatrixVal$y[2] - MatrixVal$y[1] #width of cells (h)

# Run this chunk to calculate the lambdas and the elas/sens for all years/all sites and each site-interval
######
######

#Eval SIG Growth/Surv%Metric
# for(i in 1:nrow(ModelParams_SIG)){
#   tm=ModelParams_SIG[i,]
#   tmR=ModelParams_Regional %>% filter(REGION==)
#   (tm$pg.slp*MatrixVal$y+tm$pg.int)
# }

# Change this to save the files as a new name, so they do not overwrite anything
Name <- "ALLSIG"
MODs=ModelParams_SIG
GKlist=list(nrow(MODs))
SKlist=list(nrow(MODs))
Klist=list(nrow(MODs))
Senslist=list(nrow(MODs))
Elaslist=list(nrow(MODs))
SDFlist=list(nrow(MODs))
#MODs=MODs %>% filter(!is.na(MN.RecVal_Site))
for (i in 1:nrow(MODs)){ 
  ThisModParam <- list(rec = MODs$SecRecVal[i], 
                       g.int = MODs$pg.int[i], 
                       g.slp = MODs$pg.slp[i], 
                       g.var = MODs$pg.var[i], 
                       s.int = MODs$s.int[i], 
                       s.slp = MODs$s.slp[i])
  #rec=ThisModParam$rec
  ThisMod=bigmatrix(ModelParameters = ThisModParam,MatrixVal = MatrixVal)
  MODs$Lambda[i] <- ThisMod$lam
  GKlist[[i]] <- ThisMod$Gk
  SKlist[[i]] <- ThisMod$Sk
  Klist[[i]] <- ThisMod$K
  Senslist[[i]]<- ThisMod$sens
  Elaslist[[i]] <- ThisMod$elas
  SDFlist[[i]] <- ThisMod$v.dot.w
}

MODs %>% ggplot(aes(x=round(Lambda,1),fill=REGION))+
  geom_histogram(breaks=seq(0,1.2,by=.1))+
  geom_vline(xintercept = 1)+facet_wrap("GENUS_CODE")

MODs$REPLACEMENT=round(MODs$Lambda,2)>=1

save(list=c("MODs", "GKlist","SKlist","Klist","Senslist","Elaslist","SDFlist"),
     file = "./Data/ModelData/HA_MA_AS_SIG_IMP_OUTPUTS.rdata")



#########################################################################################
plot(y,colSums(Elaslist[[65]]))

ll=lapply(X=Elaslist[which(MODs$GENUS_CODE=="POSP")],FUN = colSums)

plot(A2D(10^y),ll[[1]],type="l",xlim=c(0,30))
for(ii in 2:length(ll)){
  points(A2D(10^y),ll[[ii]],type="l")
}
MODs$Lpos=ifelse(MODs$Lambda>0.99,TRUE,FALSE)
ggplot(MODs,aes(x=log2(Lambda)))+
  geom_histogram(aes(fill=SecRecValSource),binwidth=.075)+
  geom_vline(xintercept = 0)+
  facet_wrap("REGION",scales="free_y")+
  scale_x_continuous(breaks=log2(c(.1,.33,.66,1,1.25)),labels=c(.1,.33,.66,1,1.25))



#####################################################################
#POCS
#####################################################################
#PLOT growth kernel

##Growth KERNEL
png(filename = paste0("figs/",Name,"_POCS_GrowthK.png"))
image(y, y, t(pocs$Gk)^.1, main="Pocillopora Growth Kernel",
      xlab="Colony Area T (cm^2)",ylab="Colony Area T+1 (cm^2)",xaxt="n",yaxt="n")
axis(1, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
axis(2, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
#title(GeneraNames[1], line=-1)
abline(0, 1, lty=2,lwd=3,col="black")
legend("topleft", bty = "n", lty = c(0,0,2), col = 1,
       legend= c(sprintf("recruitment value = %.2f", recvalpc),
                 sprintf("lambda = %.2f",round(pocs$lam,2)),
                 "1:1 line"))
dev.off()

##Surival KERNEL
png(filename = paste0("figs/",Name,"_POCS_SurvK.png"))
plot(y, pocs$Sk, main="Pocillopora Survivorship",type="l",
     xlab="Colony Area T (cm^2)",ylab="Prob. of Survival",ylim=c(0,1),xaxt="n")
axis(1, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
#title(GeneraNames[1], line=-1)
#abline(0, 1, lty=2)
legend("topleft", bty = "n", lty = c(0,0), col = 1,
       legend= c(sprintf("recruitment value = %.2f", recvalpc),
                 sprintf("lambda = %.2f",round(pocs$lam,2))))
dev.off()

##Growth + Survivorship KERNEL
png(filename = paste0("figs/",Name,"_POCS_PK.png"))
image(y, y, t(pocs$Pk)^.1, main="Pocillopora Growth & Survivorship Kernel",
      xlab="Colony Area T (cm^2)",ylab="Colony Area T+1 (cm^2)",xaxt="n",yaxt="n")
axis(1, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
axis(2, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
#title(GeneraNames[1], line=-1)
abline(0, 1, lty=2,lwd=3,col="black")
legend("topleft", bty = "n", lty = c(0,0,2), col = 1,
       legend= c(sprintf("recruitment value = %.2f", recvalpc), sprintf("lambda = %.2f",  
                                                                        round(pocs$lam,2)), "1:1 line"))
dev.off()

##FULL KERNEL
png(filename = paste0("figs/",Name,"_POCS_fullK.png"))
image(y, y, t(pocs$K)^.1, main="Pocillopora Full IPM Kernel",
      xlab="Colony Area T (cm^2)",ylab="Colony Area T+1 (cm^2)",xaxt="n",yaxt="n")
axis(1, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
axis(2, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
#title(GeneraNames[1], line=-1)
abline(0, 1, lty=2,lwd=3,col="black")
legend("topleft", bty = "n", lty = c(0,0,2), col = 1,
       legend= c(sprintf("recruitment value = %.2f", recvalpc), sprintf("lambda = %.2f",  
                                                                        round(pocs$lam,2)), "1:1 line"))
dev.off()
# 
# #Another way to plot full kernel
# png(filename = paste0("figs/",Name,"_POCS_alternatefullK.png"))
# image(y, y, t(pocs$K)^.1, xlab="Size (t)",ylab="Size (t+1)",col=topo.colors(100), main="Pocillopora Kernel")
# contour(y,y,t(pocs$K)^.1, add = TRUE, drawlabels = TRUE )
# legend("bottomright", bty = "n", lty = c(0,0,2), col = 1,
#        legend= c(sprintf("recruitment value = %.2f", recvalmp), sprintf("lambda = %.2f",  
#                          round(pocs$lam,2)), "1:1 line"))
# dev.off()


#Plot Elasticity
png(filename = paste0("figs/",Name,"_POCS_elasticity.png"))
#toggle having or dropping the I y[] --> y[I],pocs$elas[,]-->pocs$elas[I,]
i=2
image(y,y[],t(Elaslist[[i]])^0.2, xlab="Size (t)",ylab="Size (t+1)",main="Elasticity",xaxt="n",yaxt="n")
axis(1, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
axis(2, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000")) #[I] cuts out recruitment
abline(0,1,lwd=3,lty=2)	# plot 1:1
dev.off()

plot(y,SKlist[[6]])


#Plot Sensitivity
png(filename = paste0("figs/",Name,"_POCS_sensitivity.png"))
image(y,y,t(pocs$sens)^0.2, xlab="Size (t)",ylab="Size (t+1)", main="Pocillopora Sensitivity",xaxt="n",yaxt="n")
axis(1, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
axis(2, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
abline(0,1,lwd=3,lty=2)	# plot 1:1
dev.off()


#####################################################################
#POSP
#####################################################################
##Growth KERNEL
png(filename = paste0("figs/",Name,"_POSP_GrowthK.png"))
image(y, y, t(posp$Gk)^.1, main="Porites Growth Kernel",
      xlab="Colony Area T (cm^2)",ylab="Colony Area T+1 (cm^2)",xaxt="n",yaxt="n")
axis(1, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
axis(2, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
#title(GeneraNames[1], line=-1)
abline(0, 1, lty=2,lwd=3)
legend("topleft", bty = "n", lty = c(0,0,2), col = 1,
       legend= c(sprintf("recruitment value = %.2f", recvalpp),
                 sprintf("lambda = %.2f",round(posp$lam,2)),
                 "1:1 line"))
dev.off()

##Surival KERNEL
png(filename = paste0("figs/",Name,"_POSP_SurvK.png"))
plot(y, posp$Sk, main="Porites Survivorship",type="l",
     xlab="Colony Area T (cm^2)",ylab="Prob. of Survival",ylim=c(0,1),xaxt="n")
axis(1, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
#title(GeneraNames[1], line=-1)
#abline(0, 1, lty=2)
legend("topleft", bty = "n", lty = c(0,0), col = 1,
       legend= c(sprintf("recruitment value = %.2f", recvalpp),
                 sprintf("lambda = %.2f",round(posp$lam,2))))
dev.off()

##Growth + Survivorship KERNEL
png(filename = paste0("figs/",Name,"_POSP_PK.png"))
image(y, y, t(posp$Pk)^.1, main="Porites Growth & Survivorship Kernel",
      xlab="Colony Area T (cm^2)",ylab="Colony Area T+1 (cm^2)",xaxt="n",yaxt="n")
axis(1, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
axis(2, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
#title(GeneraNames[1], line=-1)
abline(0, 1, lty=2,lwd=3)
legend("topleft", bty = "n", lty = c(0,0,2), col = 1,
       legend= c(sprintf("recruitment value = %.2f", recvalpp),
                 sprintf("lambda = %.2f",round(posp$lam,2)),
                 "1:1 line"))
dev.off()

##FULL KERNEL
png(filename = paste0("figs/",Name,"_POSP_fullK.png"))
image(y, y, t(posp$K)^.1,main="Porites Full IPM Kernel",
      xlab="Colony Area T (cm^2)",ylab="Colony Area T+1 (cm^2)",xaxt="n",yaxt="n")
axis(1, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
axis(2, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
#title(GeneraNames[2], line=-1)
abline(0, 1, lty=2,lwd=3)
legend("topleft", bty = "n", lty = c(0,0,2), col = 1,
       legend= c(sprintf("recruitment value = %.2f", recvalpp), sprintf("lambda = %.2f",  round(posp$lam,2)), "1:1 line"))
dev.off()
# 
# #Another way to plot full kernel
# png(filename = paste0("figs/",Name,"_POSP_alternatefullK.png"))
# image(y, y, t(posp$K)^.1, xlab="Size (t)",ylab="Size (t+1)",col=topo.colors(100), main="Porites Kernel")
# contour(y,y,t(posp$K)^.1, add = TRUE, drawlabels = TRUE )
# legend("bottomright", bty = "n", lty = c(0,0,2), col = 1,
#        legend= c(sprintf("recruitment value = %.2f", recvalmp), sprintf("lambda = %.2f",  
#                          round(posp$lam,2)), "1:1 line"))
# dev.off()

#Plot Elasticity
png(filename = paste0("figs/",Name,"_POSP_elasticity.png"))
#toggle having or dropping the I y[] --> y[I],pocs$elas[,]-->pocs$elas[I,]
image(y,y[],t(posp$elas[,])^0.2, xlab="Size (t)",ylab="Size (t+1)",main="Porites Elasticity",xaxt="n",yaxt="n")
axis(1, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
axis(2, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000")) #[I] cuts out recruitment
abline(0,1,lwd=3)	# plot 1:1
dev.off()

#Plot Sensitivity
png(filename = paste0("figs/",Name,"_POSP_sensitivity.png"))
image(y,y,t(posp$sens)^0.2, xlab="Size (t)",ylab="Size (t+1)", main="Porites Sensitivity",xaxt="n",yaxt="n")
axis(1, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
axis(2, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
abline(0,1,lwd=3,lty=2)	# plot 1:1
dev.off()


#####################################################################
#MOSP
#####################################################################
##Growth KERNEL
png(filename = paste0("figs/",Name,"_MOSP_GrowthK.png"))
image(y, y, t(mosp$Gk)^.1, main="Montipora Growth Kernel",
      xlab="Colony Area T (cm^2)",ylab="Colony Area T+1 (cm^2)",xaxt="n",yaxt="n")
axis(1, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
axis(2, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
#title(GeneraNames[1], line=-1)
abline(0, 1, lty=2,lwd=3)
legend("topleft", bty = "n", lty = c(0,0,2), col = 1,
       legend= c(sprintf("recruitment value = %.2f", recvalmp),
                 sprintf("lambda = %.2f",round(mosp$lam,2)),
                 "1:1 line"))
dev.off()

##Surival KERNEL
png(filename = paste0("figs/",Name,"_MOSP_SurvK.png"))
plot(y, mosp$Sk, main="Montipora Survivorship",type="l",
     xlab="Colony Area T (cm^2)",ylab="Prob. of Survival",ylim=c(0,1),xaxt="n")
axis(1, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
#title(GeneraNames[1], line=-1)
#abline(0, 1, lty=2)
legend("topleft", bty = "n", lty = c(0,0), col = 1,
       legend= c(sprintf("recruitment value = %.2f", recvalmp),
                 sprintf("lambda = %.2f",round(mosp$lam,2))))
dev.off()

##Growth + Survivorship KERNEL
png(filename = paste0("figs/",Name,"_MOSP_PK.png"))
image(y, y, t(mosp$Pk)^.1, main="Montipora Growth & Survivorship Kernel",
      xlab="Colony Area T (cm^2)",ylab="Colony Area T+1 (cm^2)",xaxt="n",yaxt="n")
axis(1, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
axis(2, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
#title(GeneraNames[1], line=-1)
abline(0, 1, lty=2,lwd=3)
legend("topleft", bty = "n", lty = c(0,0,2), col = 1,
       legend= c(sprintf("recruitment value = %.2f", recvalmp),
                 sprintf("lambda = %.2f",round(mosp$lam,2)),
                 "1:1 line"))
dev.off()

##FULL KERNEL
png(filename = paste0("figs/",Name,"_MOSP_fullK.png"))
image(y, y, t(mosp$K)^.1,main="Montipora Full IPM Kernel",
      xlab="Colony Area T (cm^2)",ylab="Colony Area T+1 (cm^2)",xaxt="n",yaxt="n")
axis(1, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
axis(2, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
#title(GeneraNames[2], line=-1)
abline(0, 1, lty=2,lwd=3)
legend("topleft", bty = "n", lty = c(0,0,2), col = 1,
       legend= c(sprintf("recruitment value = %.2f", recvalmp), sprintf("lambda = %.2f",  
                                                                        round(mosp$lam,2)), "1:1 line"))
dev.off()
# 
# #Another way to plot full kernel
# png(filename = paste0("figs/",Name,"_MOSP_alternatefullK.png"))
# image(y, y, t(mosp$K)^.1, xlab="Size (t)",ylab="Size (t+1)",col=topo.colors(100), main="Montipora Kernel")
# contour(y,y,t(mosp$K)^.1, add = TRUE, drawlabels = TRUE )
# legend("bottomright", bty = "n", lty = c(0,0,2), col = 1,
#        legend= c(sprintf("recruitment value = %.2f", recvalmp), sprintf("lambda = %.2f",  
#                          round(mosp$lam,2)), "1:1 line"))
# dev.off()

#Plot Elasticity
png(filename = paste0("figs/",Name,"_MOSP_elasticity.png"))
#toggle having or dropping the I y[] --> y[I],pocs$elas[,]-->pocs$elas[I,]
image(y,y[],t(mosp$elas[,])^0.2, xlab="Size (t)",ylab="Size (t+1)",main="Montipora Elasticity",xaxt="n",yaxt="n")
axis(1, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
axis(2, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000")) #[I] cuts out recruitment
abline(0,1,lwd=3)	# plot 1:1
dev.off()

#Plot Sensitivity
png(filename = paste0("figs/",Name,"_MOSP_sensitivity.png"))
image(y,y,t(mosp$sens)^0.2, xlab="Size (t)",ylab="Size (t+1)", main="Montipora Sensitivity",xaxt="n",yaxt="n")
axis(1, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
axis(2, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
abline(0,1,lwd=3,lty=2)	# plot 1:1
dev.off()


##All Three Survival KERNEL
png(filename = paste0("figs/",Name,"_POCS_SurvK.png"))
plot(y, posp$Sk, main="Three Taxa Comparasion: Survivorship",type="l",
     xlab="Colony Area T (cm^2)",ylab="Prob. of Survival",ylim=c(0,1),xaxt="n",col="darkgreen")
points(y, mosp$Sk,type="l",col="darkblue")
points(y, pocs$Sk,type="l",col="darkred")
axis(1, at = seq(0, 5, 1), labels = c("1","10","100","1,000","10,000","100,000"))
#title(GeneraNames[1], line=-1)
#abline(0, 1, lty=2)
legend("topleft", bty = "n", lty = c(0,0), col = 1,
       legend= c(sprintf("recruitment value = %.2f", recvalpc),
                 sprintf("lambda = %.2f",round(pocs$lam,2))))
dev.off()



#Plot stable size distribution
#png(filename = "figs/AllInt_MOSP_stablesizedist.png")
plot(y, mosp$w,xlab="Size",type="l",main="MOSP Stable size distribution")
subset(ColonyLevel$Genus_Code=="MOSP")#dev.off()

#Plot left eigenvector (reproductive value)
#plot(y,reprodvalue,xlab="Size",type="l",main="Reproductive values")



png("Figures/VRmodel_Params.png", height = 30*nrow(data), width = 90*ncol(data))
grid.table(data)
dev.off()

```

#Perturbation Analysis
```{r PerturbationAnalysis, }
#Assessing the sensitivity of each component model parameter individually, by adjusting each coefficient 
# by +/-1%, and see how it affect the population growth rate (lambda)

#change parameter list to dataframe
df_MP_POSP <- ldply(MP_POSP, data.frame)
modparam_POSP <- spread(df_MP_POSP, key = .id, value = X..i.. )
df_MP_POCS <- ldply(MP_POCS, data.frame)
modparam_POCS <- spread(df_MP_POCS, key = .id, value = X..i.. )
df_MP_MOSP <- ldply(MP_MOSP, data.frame)
modparam_MOSP <- spread(df_MP_MOSP, key = .id, value = X..i.. )



ChangeMe <- modparam_POSP #!!!!!

nparams <- length(ChangeMe) #total number of component model parameters
lam.sn <- matrix(ncol = nparams, nrow = 2, data = NA)  #matrix where every row
all.params.orig <- ChangeMe

for (p in 1:nparams) { #tweak each parameter in turn
  for (j in 1:2) {     #first down by 1% and then 1% up
    all.params <- all.params.orig
    all.params[p] <- all.params[p]*(1+2*(j-1.5)/100)
    BG <- bigmatrix(all.params)
    lam.sn[j,p] <-BG$lam
  }
}
lam.sn #top row is decr by 1%, bottom row is incr

#POSP
normPOSP <- matrix(posp$lam, nrow = 2, ncol = 6)
(round((normPOSP -lam.sn)/posp$lam, 3)) *100

#POCS
normPOCS <- matrix(pocs$lam, nrow = 2, ncol = 6)
(round((normPOCS - lam.sn)/pocs$lam, 3)) *100

#MOSP
normMOSP <- matrix(mosp$lam, nrow = 2, ncol = 6)
(round((normMOSP -lam.sn)/posp$lam, 3)) *100


# perturb manually for -0.1 and then +0.1:
all.params <- all.params.orig
all.params[10] <- -0.01
growth2 <- delta_size * outer(y, y, growth_function, g.int=POSP$g.int[1],   #growth kernel
                              g.slp=POSP$g.slp[1], g.var=POSP$g.var[1])
survival2 <- survival_function(y, s.int=POSP$s.int[1], s.slp=POSP$s.slp[1])     #survival kernel
P2 <- Growth_kernel    #placeholder
z <- 1:n
P2[,z] <- Growth_kernel[,z]*Survival_kernel[z]    #growth/survival kernel
# P2 <- growth2 * survival2
# for (z in 1:n) {
#   P2[,z] <- growth2[,z]*survival2[z]
#   }
Repro2 <- delta_size * outer(y, y, reproduction_function_josh)    #reproduction kernel
K2 <- P2 + Repro2
lam_new <- Re(eigen(K2)$values[1])
lam.sn[21] <- lam_new

```




The next chunk is for the plots Tom made
```{r TomsPlots,}
Name=Scenario
load(file = sprintf("data/%s_all_Lambdas.rdata",Name))

AllD=MODs#merge(data,Lambdas)
#AllD$AY=c(T,T,T,rep(F,nrow(AllD)-3))
#AllD$Island=substr(AllD$Sites,1,3)
#names(AllD)[which(names(AllD)=="AIC")]=c("g.AIC","s.AIC")
AllD$Sites=factor(AllD$Sites,levels=
                    sort(unique(AllD$Sites))[c(1,9,16,2,3,14,15,10:13,4:8)])
AllD$Years=factor(AllD$Years,levels=
                    sort(unique(AllD$Years)))
AllD$EndYear=2000+as.numeric(substr(AllD$Years,4,5))
AllD$EndYear[AllD$Years=="All Years"&AllD$Genus_Code=="POSP"]=max(AllD$EndYear[AllD$Genus_Code=="POSP"],na.rm=T)+.5
AllD$EndYear[AllD$Years=="All Years"&AllD$Genus_Code=="POCS"]=max(AllD$EndYear[AllD$Genus_Code=="POCS"],na.rm=T)+.5
AllD$EndYear[AllD$Years=="All Years"&AllD$Genus_Code=="MOSP"]=max(AllD$EndYear[AllD$Genus_Code=="MOSP"],na.rm=T)+.5
AllD$StartYear=2000+as.numeric(substr(AllD$Years,1,2))
AllD$StartYear[AllD$Years=="All Years"&AllD$Genus_Code=="POSP"]=min(AllD$StartYear[AllD$Genus_Code=="POSP"],na.rm=T)-.5
AllD$StartYear[AllD$Years=="All Years"&AllD$Genus_Code=="POCS"]=min(AllD$StartYear[AllD$Genus_Code=="POCS"],na.rm=T)-.5
AllD$StartYear[AllD$Years=="All Years"&AllD$Genus_Code=="MOSP"]=min(AllD$StartYear[AllD$Genus_Code=="MOSP"],na.rm=T)-.5
AllD$Island=factor(AllD$Island,levels=c("All","HAW","MAI","OAH","FFS","KUR"))
GeneraLU=c("Porites sp.","Montipora sp.","Pocillopora sp.")
names(GeneraLU)=c("POSP","MOSP","POCS")
AllD$Genus=GeneraLU[AllD$Genus_Code]
AllD$Genus=factor(AllD$Genus,levels=GeneraLU)
#subset(AllD,Island%in%c("All","HAW","MAI","OAH"))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(length(unique(AllD$Sites)))
tobelog2Breaks=c(.5,.75,1,1.5,2,3)#seq(.25,3.5,by=.25)
log2Breaks=log2(tobelog2Breaks)
LambdaBySite=ggplot(subset(AllD,Sites!="All Sites"),aes(x=EndYear,y=log2(Lambdas),size=GModel.n,fill=log2(Lambdas)))+
  geom_segment(aes(xend=StartYear,yend=log2(Lambdas)),alpha=.5,size=1.5)+
  geom_path(color="gray50",lty=2,size=1)+
  geom_point(aes(shape=Island),color="black")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 2014.75,color="orange")+
  geom_vline(xintercept = 2015.75,color="red")+
  geom_vline(xintercept = 2019.75,color="orange")+
  facet_wrap(Genus~Sites,ncol=7)+
  scale_color_manual(name="Site",values = mycolors)+
  scale_fill_gradient2(name="Lambda",midpoint = 0,breaks = log2Breaks,labels = tobelog2Breaks)+
  scale_size_binned(name="N Colonies",trans = "log10")+
  scale_y_continuous(breaks = log2Breaks,labels = tobelog2Breaks,name="Lambda - Const. Region Recruitment")+
  scale_shape_manual(name="Island",values = c(20:25))+
  theme_bw()+theme(legend.position = "bottom")
LambdaBySite
ggsave(filename = "./figs/LambdasBySitePlots.jpg",plot = LambdaBySite,width=16,height=9)

LambdaByIsland=ggplot(AllD,aes(x=EndYear,y=log2(Lambdas),size=GModel.n,fill=log2(Lambdas)))+
  geom_segment(aes(xend=StartYear,yend=log2(Lambdas)),alpha=.5,size=1.5)+
  geom_path(aes(color=Sites),size=1)+
  geom_point(aes(shape=Island),color="black")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 2014.75,color="orange")+
  geom_vline(xintercept = 2015.75,color="red")+
  geom_vline(xintercept = 2019.75,color="orange")+
  facet_grid(Island~Genus)+
  scale_color_manual(name="Site",values = mycolors)+
  scale_fill_gradient2(name="Lambda",midpoint = 0,breaks = log2Breaks,labels = tobelog2Breaks)+
  scale_size_binned(name="N Colonies",trans = "log10")+
  scale_y_continuous(breaks = log2Breaks,labels = tobelog2Breaks,name="Lambda - Const. Region Recruitment")+
  scale_shape_manual(name="Island",values = c(20:25))+
  theme_bw()
LambdaByIsland
ggsave(filename = "./figs/LambdaByIslandPlots.jpg",plot = LambdaByIsland,width=8,height=8)


LambdaAll=ggplot(AllD,aes(x=EndYear,y=log2(Lambdas),size=GModel.n,fill=log2(Lambdas)))+
  geom_segment(aes(xend=StartYear,yend=log2(Lambdas)),alpha=.5,size=1.5)+
  #geom_path(aes(color=Sites),size=1)+
  geom_point(aes(shape=Island),color="black")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 2014.75,color="orange")+
  geom_vline(xintercept = 2015.75,color="red")+
  geom_vline(xintercept = 2019.75,color="orange")+
  facet_grid(.~Genus)+
  scale_color_manual(name="Site",values = mycolors)+
  scale_fill_gradient2(name="Lambda",midpoint = 0,breaks = log2Breaks,labels = tobelog2Breaks)+
  scale_size_binned(name="N Colonies",trans = "log10")+
  scale_y_continuous(breaks = log2Breaks,labels = tobelog2Breaks,name="Lambda - Const. Region Recruitment")+
  scale_shape_manual(name="Island",values = c(20:25))+
  theme_bw()
LambdaAll
ggsave(filename = "./figs/LambdasAllPlots.jpg",plot = LambdaAll,width=12,height=6)

names(AllD)
demo=AllD[,c("Island","Sites","Years","StartYear","EndYear","Genus","GModel.n","g.int","g.slp","g.var","s.int","s.slp","Lambdas")]
demo_num=demo[,8:12]

library(vegan)
library(ggfortify)
pca_res <- prcomp(demo_num, scale = TRUE)
demo=cbind(demo,pca_res$x)
sc=4
arrow_pca=as.data.frame(pca_res$rotation)
arrow_pca$name=row.names(arrow_pca)
hc=hclust(dist(demo_num))

h2k=function(hc,hrange){
  h2k_out=data.frame(h=hrange,k=NA)
  for(i in 1:length(hrange)){
    h2k_out$k[i]=length(unique(cutree(hc,h=hrange[i])))
  }
  return(h2k_out)
}
hr=seq(0,40,by=.1)
k=h2k(hc,hr)
sort(table(k$k),decreasing = T)

ID=cutree(hc,k = 4)
demo$ClustID=factor(ID)
plot(hc)
abline(h = 4.5,col="red")
demo$AY=demo$Years=="All Years"
ggplot(demo,aes(x=PC2,y=PC3))+
  geom_segment(aes(x = 0,xend=sc*PC2,y=0,yend=sc*PC3),data=arrow_pca,color="gray50")+
  geom_text(aes(x=sc*PC2,y=sc*PC3,label=name),data=arrow_pca,color="gray50")+
  geom_point(data=subset(demo,AY==T),aes(color=Genus),size=5,pch=1)+
  geom_text(aes(shape=ClustID,color=Genus,label=ClustID,size=Lambdas))+
  scale_shape_manual(values=1:15)+
  scale_color_discrete()+
  theme_bw()+facet_grid("Genus")

plot(log10(demo$Lambdas)~demo$ClustID)
abline(h=0)

library(scales)
ggplot(AllD,aes(group=Genus_Code))+
  geom_hline(yintercept = c(14.75),size=3,color="gold")+
  geom_hline(yintercept = c(15.75),size=3,color="red")+
  geom_hline(yintercept = c(19.75),size=3,color="orange")+
  geom_linerange(aes(ymin=StartYear,ymax=EndYear,x=Sites,color=log2(Lambdas),size=log2(Lambdas)),
                 position=position_dodge(width=.5),alpha=1)+#,data=subset(AllD,Genus_Code=="POCS"))+
  geom_point(aes(shape=Genus_Code,y=StartYear,x=Sites),fill="black",size=2,position=position_dodge(width=.5),color="black",alpha=1)+#,data=subset(AllD,Genus_Code=="POCS"))+
  geom_point(aes(shape=Genus_Code,y=EndYear,x=Sites,size=log2(Lambdas),fill=log2(Lambdas)),position=position_dodge(width=.5),color="black",alpha=1)+
  geom_vline(xintercept = c(1.5,2.5,3.5,5.5,7.5,11.5))+
  scale_fill_distiller(type = "div",name="Lambda",palette = "Spectral",direction = 1,limits = c(-1,1))+
  #scale_fill_gradient2(name="Lambda",mid = "yellow",high="blue",low="red")+
  scale_size_area(name="Lambda")+
  scale_color_distiller(type = "div",name="Lambda",palette = "Spectral",direction = 1,limits = c(-1,1))+
  #scale_color_gradient2(name="Lambda",mid = "yellow",high="blue",low="red")+
  scale_shape_manual(name="Genus",values = c(22,21,23))+
  #  scale_size_area(limits=log2(c(.5,2)),breaks=log2(seq(.5,2,length.out = 5)),labels=seq(.5,2,length.out = 5),oob=squish)+
  coord_equal()+
  coord_flip()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,size = 10,hjust = 0),panel.background=element_rect(fill = "white"),panel.grid = element_blank())

ggplot(AllD,aes(group=Genus_Code))+
  geom_hline(yintercept = c(14.75),size=3,color="gold")+
  geom_hline(yintercept = c(15.75),size=3,color="red")+
  geom_hline(yintercept = c(19.75),size=3,color="orange")+
  geom_linerange(aes(ymin=StartYear,ymax=EndYear,x=Sites,color=log2(Lambdas),size=log2(Lambdas)),
                 position=position_dodge(width=.5),alpha=1)+#,data=subset(AllD,Genus_Code=="POCS"))+
  geom_point(aes(shape=Genus_Code,y=StartYear,x=Sites),fill="black",size=2,position=position_dodge(width=.5),color="black",alpha=1)+#,data=subset(AllD,Genus_Code=="POCS"))+
  geom_point(aes(shape=Genus_Code,y=EndYear,x=Sites,size=log2(Lambdas),fill=log2(Lambdas)),position=position_dodge(width=.5),color="black",alpha=1)+
  geom_vline(xintercept = c(1.5,2.5,3.5,5.5,7.5,11.5))+
  scale_fill_distiller(type = "div",name="Lambda",palette = "Spectral",direction = 1,limits = c(-1,1))+
  #scale_fill_gradient2(name="Lambda",mid = "yellow",high="blue",low="red")+
  scale_size_area(name="Lambda")+
  scale_color_distiller(type = "div",name="Lambda",palette = "Spectral",direction = 1,limits = c(-1,1))+
  #scale_color_gradient2(name="Lambda",mid = "yellow",high="blue",low="red")+
  scale_shape_manual(name="Genus",values = c(22,21,23))+
  #  scale_size_area(limits=log2(c(.5,2)),breaks=log2(seq(.5,2,length.out = 5)),labels=seq(.5,2,length.out = 5),oob=squish)+
  coord_equal()+
  coord_flip()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,size = 10,hjust = 0),panel.background=element_rect(fill = "white"),panel.grid = element_blank())



ggplot(subset(AllD,Sites=="MAI_SIO_K02"),aes(x=log2(Lambdas)))+
  geom_histogram(aes(fill=R.sqd,color=pR2))+
  geom_vline(aes(xintercept = log2(Lambdas)),color="red",data=subset(AllD,Sites=="All Sites"))+
  facet_grid(.~Genus_Code)+
  theme_bw()



```
Ignore the stuff below!!
  
  
  ```{r LoadData_Subsetting,}
# Loading in the data from analysis.Rmd
load("data/Growth_intmod_params.rdata")
load("data/Survival_intmod_params.rdata")
# load("data/Recruitment_sizeIndependent.rdata")
# load("data/Recruitment_sizeDependent.rdata")
load("data/SiteInterval_Gmodfits.rdata")
load("data/SiteInterval_Smodfits.rdata")

# Separating the growth data frame by genera #
POCS_G_1 <- dt[dt$Genus_Code == "POCS",][1,]
POCS_G_1 <- cbind(POCS_G_1, data.frame(g.var=0.230512176))
a        <- GrowthDataFrame[GrowthDataFrame$Genus_Code == "POCS",]
POCS_G_2 <- a[which(a$g.int >0),1:6]
POCS_G   <- rbind(POCS_G_1, POCS_G_2) 

POSP_G_1   <- dt[dt$Genus_Code == "POSP",][1,]
POSP_G_1 <- cbind(POSP_G_1, data.frame(g.var=0.372585621))
b        <- GrowthDataFrame[GrowthDataFrame$Genus_Code == "POSP",]
POSP_G_2 <- b[which(b$g.int >0),1:6]
POSP_G   <- rbind(POSP_G_1, POSP_G_2) 

MOSP_G_1   <- dt[dt$Genus_Code == "MOSP",][1,]
MOSP_G_1 <- cbind(MOSP_G_1, data.frame(g.var=1.102913776))
c        <- GrowthDataFrame[GrowthDataFrame$Genus_Code == "MOSP",] 
MOSP_G_2 <- c[which(c$g.int >0),1:6]
MOSP_G   <- rbind(MOSP_G_1, MOSP_G_2)

##############
#SURVIVAL
# subsetting the survival data frame by genera 
POCS_S_1 <- dS[dS$Genus_Code == "POCS",][1,]
a        <- SurvivalDataFrame[SurvivalDataFrame$Genus_Code == "POCS",]
POCS_S_2 <- a[which(a$s.int >0),1:5]
POCS_S   <- rbind(POCS_S_1,POCS_S_2)

POSP_S_1 <- dS[dS$Genus_Code == "POSP",][1,]
b        <- SurvivalDataFrame[SurvivalDataFrame$Genus_Code == "POSP",]
POSP_S_2 <- b[which(b$s.int >0),1:5]
POSP_S   <- rbind(POSP_S_1,POSP_S_2)

MOSP_S_1 <- dS[dS$Genus_Code == "MOSP",][1,]
c        <- SurvivalDataFrame[SurvivalDataFrame$Genus_Code == "MOSP",] 
MOSP_S_2 <- c[which(c$s.int >0),1:5]
MOSP_S   <- rbind(MOSP_S_1,MOSP_S_2)


# ##############
# #RECRUITMENT
# POCS_R <- subset(dr, Genus_Code == "POCS")
# POSP_R <- subset(dr, Genus_Code == "POSP")
# MOSP_R <- subset(dr, Genus_Code == "MOSP")
# 
# POCS_Rdep <- subset(dr_size, Genus_Code == "POCS")
# POSP_Rdep <- subset(dr_size, Genus_Code == "POSP")
# MOSP_Rdep <- subset(dr_size, Genus_Code == "MOSP")
# 
# POCS_R <- subset(POCS_R,Years == "All_Years" )
# totrecruits_POCS <- sum(POCS_R$Tots_recs_peryear, na.rm = T)
# mean(POCS_R$Tots_recs_peryear, na.rm = T)
# 
# POSP_R <- subset(POSP_R,Years == "All_Years" )
# totrecruits_POSP <- sum(POSP_R$Tots_recs_peryear, na.rm = T)
# mean(POSP_R$Tots_recs_peryear, na.rm = T)
# 
# MOSP_R <- subset(MOSP_R,Years == "All_Years" )
# totrecruits_MOSP <- sum(MOSP_R$Tots_recs_peryear, na.rm = T)
# mean(MOSP_R$Tots_recs_peryear, na.rm = T)


```

