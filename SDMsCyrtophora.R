###Code for SDMs in Segura et.al_Cyrtophora invasion in America_BIOL.INV####
rm(list=ls())
gc()

setwd("C:/Users/Pichi/Dropbox/PANDEMIC PRIORITIES/ReviewCyrtophoraCitricola/")

#Packages
PKGs <- c("raster", "ENMeval", "dplyr","rJava", "rgdal","jsonlite", "maptools", 
          "sp", "spThin", "usdm", "sf", "ggplot2")

#Load libraries
sapply(PKGs, require, character.only = TRUE)

#Load data
#Climatic Variables Resolution 2.5 arcmins

##WORLDCLIM
WC2.5<-raster::getData('worldclim', var='bio', res=2.5)

###MAX WIND VELOCITY
wind<-raster('wc2.1_2.5m_wind/windmax.asc')
wind<-resample(wind, WC2.5[[1]])

###STACK ALL VARIABLES  
AllVars<-stack(WC2.5[[10]],WC2.5[[16]],WC2.5[[17]],WC2.5[[2]],WC2.5[[15]],wind) ##selected
plot(AllVars[[1]])

#Occurrences
amer<-read.csv("Datos/Occs/todosamerica.csv")
points(amer[,3:2], col="red", pch=19)

medit<-read.csv("Datos/Occs/todosmediterraneo.csv")
points(medit[,3:2], col="darkgreen", pch=19)

suraf<-read.csv("Datos/Occs/todossurafrica.csv")
points(suraf[,3:2], col="blue", pch=19)

##Projection areas
eAM<-c(-128,-30,-47,52) 
AME<-raster::crop(AllVars,eAM)

eMED<-c(-12,40,20,50)
MED<-crop(AllVars,eMED)

eSUD<-c(10,41,-40,-15)
SUD<-crop(AllVars,eSUD)

#50km buffers to define calibration areas
BuffAmer<-readOGR("Datos", layer = "BuffAM50km")
BuffMedit<-readOGR("Datos", layer = "BuffMED50km")
BuffSuraf<-readOGR("Datos", layer = "BuffSUD50km")

###Calibration areas
bg_AME<-mask(AME,BuffAmer)
bg_AME<-crop(bg_AME,BuffAmer)
plot(bg_AME[[1]])

bg_MED<-mask(MED,BuffMedit)
bg_MED<-crop(bg_MED,BuffMedit)
plot(bg_MED[[1]])

bg_SUD<-mask(SUD,BuffSuraf)
bg_SUD<-crop(bg_SUD,BuffSuraf)
plot(bg_SUD[[1]])

###THINNING
am_OCCS_5km<-as.data.frame(thin.algorithm(rec.df.orig = amer,thin.par = 5,reps = 1))
names(am_OCCS_5km)<-names(amer)

med_OCCS_5km<-as.data.frame(thin.algorithm(rec.df.orig = medit,thin.par = 5,reps = 1))
names(med_OCCS_5km)<-names(medit)

sud_OCCS_5km<-as.data.frame(thin.algorithm(rec.df.orig = suraf,thin.par = 5,reps = 1))
names(sud_OCCS_5km)<-names(suraf)

### BACKGROUND RANDOM POINTS (10000per region)
bgpoints_AME <- dismo::randomPoints(bg_AME[[1]], n = 10000) %>% as.data.frame()
names(bgpoints_AME)<-c("Longitude", "Latitude")

bgpoints_MED <- dismo::randomPoints(bg_MED[[1]], n = 10000) %>% as.data.frame()
names(bgpoints_MED)<-c("Longitude", "Latitude")

bgpoints_SUD <- dismo::randomPoints(bg_SUD[[1]], n = 10000) %>% as.data.frame()
names(bgpoints_SUD)<-c("Longitude", "Latitude")
points(bgpoints_SUD)


###################################################################################
#################                CANDIDATE MODELS               ###################
###################################################################################

####AMERICA
AME_mods <- ENMevaluate(updateProgress = TRUE,occs = am_OCCS_5km[,3:2], envs = bg_AME, bg = bgpoints_AME, 
                        algorithm = 'maxent.jar', partitions = 'checkerboard2',parallel = TRUE, 
                        progbar = TRUE,tune.args = list(fc = c("L","H","Q","LH","LQ","QH","LQH"), 
                        rm = seq(from=1, to=4, by=0.5)), overlap = FALSE)

###RESULTS
res_AME <- eval.results(AME_mods)

###OPTIMAL AUC & 10p OR
opt.seqAME <- res_AME %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seqAME

###CHOSEN MODEL
OptMod_AME <- eval.models(AME_mods)[[opt.seqAME$tune.args]]

###CHECK RESPONSES
plot(OptMod_AME, type = "cloglog")

dev.off()
palette<-colorRampPalette(c("dodgerblue4","deepskyblue4","gold","red3"))(255)

plot(eval.predictions(AME_mods)[[opt.seqAME$tune.args]], ylim = c(-20,30), xlim = c(-92,-42), 
     legend = FALSE, col=palette,main = paste('Prediction',as.character(opt.seqAME$tune.args)))

##PROJECTION IN THE WHOLE EXTENT OF AMERICA
AME_opt<-dismo::predict(OptMod_AME, AME, args=c("outputformat=logistic"))
plot(AME_opt,col=palette)
writeRaster(AME_opt, "NEW_Models_RASTERS/AMEbg_to_AME_opt.asc")

###########################################################################################

###########################################################################################

####MEDITERRANEAN
MED_mods <- ENMevaluate(occs = med_OCCS_5km[,3:2],envs = bg_MED, bg = bgpoints_MED, 
                         algorithm = 'maxent.jar', method = 'checkerboard2',parallel = TRUE, 
                         progbar = TRUE,tune.args = list(fc = c("L","H","Q","LH","LQ","QH","LQH"), 
                         rm = seq(from=1, to=4, by=0.5)), overlap = FALSE)

###RESULTS
res_MED <- eval.results(MED_mods)

###OPTIMAL AUC & 10p OR
opt.seqMED <- res_MED %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))

opt.seqMED

###CHOSEN MODEL
OptMod_MED <- eval.models(MED_mods)[[opt.seqMED$tune.args]]

###CHECK RESPONSES
plot(OptMod_MED, type = "cloglog")

dev.off()
plot(eval.predictions(MED_mods)[[opt.seqMED$tune.args]], ylim = c(-10,70), xlim = c(-50,100), 
     legend = FALSE, col=palette,main = paste('Prediction',as.character(opt.seqMED$tune.args)))

MED_to_MED_opt<-dismo::predict(OptMod_MED, MED, args=c("outputformat=logistic"))
plot(MED_to_MED_opt,col=palette)
writeRaster(MED_to_MED_opt, "NEW_Models_RASTERS/MEDbg_to_MED_opt.asc", overwrite=TRUE)

##PROJECTION IN THE WHOLE EXTENT OF AMERICA
MED_to_AME_opt<-dismo::predict(OptMod_MED, AME, args=c("outputformat=logistic"))
plot(MED_to_AME_opt,col=palette)
writeRaster(MED_to_AME_opt, "NEW_Models_RASTERS/MEDbg_to_AME_opt.asc")

###########################################################################################

###########################################################################################

####SUDAFRICA
SUD_mods <- ENMevaluate(occs = sud_OCCS_5km[,3:2],envs = SUD, bg = bgpoints_SUD, 
                        algorithm = 'maxent.jar', method = 'checkerboard2',parallel = TRUE, 
                        progbar = TRUE,tune.args = list(fc = c("L","H","Q","LH","LQ","QH","LQH"), 
                                                        rm = seq(from=1, to=4, by=0.5), overlap = FALSE))
###RESULTS
res_SUD <- eval.results(SUD_mods)

###OPTIMAL AUC & 10p OR
opt.seqSUD <- res_SUD %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))

opt.seqSUD

###CHOSEN MODEL
OptMod_SUD <- eval.models(SUD_mods)[[opt.seqSUD$tune.args]]

###CHECK RESPONSES
plot(OptMod_SUD, type = "cloglog")

dev.off()
plot(eval.predictions(SUD_mods)[[opt.seqSUD$tune.args]], ylim = c(-40,-15), xlim = c(10,41), 
     legend = FALSE, col=palette, main = paste('Prediction',as.character(opt.seqSUD$tune.args)))

SUD_to_SUD_opt<-dismo::predict(OptMod_SUD, SUD, args=c("outputformat=logistic"))
plot(SUD_to_SUD_opt, col=palette)
writeRaster(SUD_to_SUD_opt, "NEW_Models_RASTERS/SUDbg_to_SUD_opt.asc")

####PROJECTION TO WHOLE AMERICA
SUD_to_AME_opt<-dismo::predict(OptMod_SUD, AME2, args=c("outputformat=logistic"))
plot(SUD_to_AME_opt,col=palette)
writeRaster(SUD_to_AME_opt, "NEW_Models_RASTERS/SUDbg_to_AME_opt.asc")

########################################################################################

########################################################################################

####PLOTTING GENERAL RESULTS
par(mfrow=c(2,1))
plot(eval.predictions(AME_mods)[[opt.seqAME$tune.args]], ylim = c(-20,30), xlim = c(-92,-42), 
     legend = FALSE, col=palette,main = paste('AMER (AMER POINTS)\n',as.character(opt.seqAME$tune.args)))
points(amer[,3:2], col="black",cex=0.75, pch=19)

plot(eval.predictions(MED_mods)[[opt.seqMED$tune.args]], ylim = c(-10,70), xlim = c(-50,100), 
     legend = FALSE, col=palette, main = paste('MED (MED POINTS)\n',as.character(opt.seqMED$tune.args)))
points(medit[,3:2], col="black",cex=0.75, pch=19)

plot(eval.predictions(SUD_mods)[[opt.seqSUD$tune.args]], ylim = c(-40,-15), xlim = c(10,41), 
     legend = FALSE, col=palette, main = paste('SUDAF(SUDAF POINTS)\n',as.character(opt.seqSUD$tune.args)))
points(suraf[,3:2], col="black",cex=0.75, pch=19)

plot(AME_opt,col=palette, main = paste('AMER (AMER POINTS)\n',as.character(opt.seqAME$tune.args)))
points(amer[,3:2], col="black",cex=0.75, pch=19)

plot(MED_to_AME_opt,col=palette, main = paste('AMER (MED POINTS)\n',as.character(opt.seqMED$tune.args)))
points(amer[,3:2], col="black",cex=0.75, pch=21)

plot(SUD_to_AME_opt,col=palette, main = paste('AMER (SUD POINTS)\n',as.character(opt.seqSUD$tune.args)))
points(amer[,3:2], col="black",cex=0.75, pch=21)


###DENSITY PLOTS SUITABILITY ASSIGNED TO ALIEN RECORDS
SuitSUD_AMERpoints<-as.data.frame(raster::extract(SUD_to_AME_opt,am_OCCS_5km[,3:2]))
Obs<-rep("SUD",length(SuitSUD_AMERpoints))
DF_SUD<-as.data.frame(cbind(Obs,am_OCCS_5km[,3:2],SuitSUD_AMERpoints))
names(DF_SUD)<-c("Obs","Longitude","Latitude","Suitability")
write.csv(DF_SUD,"SuitPoints_fromSUD_THIN.csv")

SuitMED_AMERpoints<-as.data.frame(raster::extract(MED_to_AME_opt,am_OCCS_5km[,3:2]))
Obs2<-rep("MED",length(SuitSUD_AMERpoints))
DF_MED<-as.data.frame(cbind(Obs2,am_OCCS_5km[,3:2], SuitMED_AMERpoints))
names(DF_MED)<-c("Obs","Longitude","Latitude","Suitability")
write.csv(DF_MED,"SuitPoints_fromMED_THIN.csv")

DFsuitVals<-rbind(DF_SUD,DF_MED)
write.csv(DFsuitVals,"DFsuitValsTHINoccsAMER.csv")

library(ggplot2)

dens_SUD<-ggplot(DF_SUD, aes(x = `Suitability`, y = Obs, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "Suitability", option = "B") +
  theme_classic()


dens_MED<-ggplot(DF_MED, aes(x = `Suitability`, y = Obs, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "Suitability", option = "B") +
  theme_classic()

library(ggpubr)
ggarrange(dens_MED,dens_SUD,ncol = 1, nrow = 2)

###Data for heatmap of points
normSUIT_MED<-(SuitMED_AMERpoints-min(SuitMED_AMERpoints))/(max(SuitMED_AMERpoints)-min(SuitMED_AMERpoints))
normSUIT_SUD<-(SuitSUD_AMERpoints-min(SuitSUD_AMERpoints))/(max(SuitSUD_AMERpoints)-min(SuitSUD_AMERpoints))

Pts_Suit<-cbind(amer[,3:2],Vals)

write.csv(Pts_Suit,"Pts_Suit.csv")

###BOYCE INDEX
boyce_MED<-ecospat.boyce (MED_to_AME_opt, amer[,3:2], nclass = 0, window.w = "default", res = 100,
                          PEplot = TRUE)$Spearman.cor    

boyce_SUD<-ecospat.boyce (SUD_to_AME_opt, amer[,3:2], nclass = 0, window.w = "default", res = 100,
                          PEplot = TRUE)$Spearman.cor    

###############################################################################################