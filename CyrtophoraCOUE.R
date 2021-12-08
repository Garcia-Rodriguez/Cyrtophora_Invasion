###ECOSPAT CYRTOPHORA
library(ecospat)

####VARIABLES
##BIOCLIM
WC2.5<-raster::getData('worldclim', var='bio', res=2.5)

###MAX WIND VELOCITY
wind<-raster('wc2.1_2.5m_wind/windmax.asc')
wind<-resample(wind, WC2.5[[1]])

###STACK VARIABLES  
SelVars<-stack(WC2.5[[10]],WC2.5[[16]],WC2.5[[17]],WC2.5[[2]],WC2.5[[15]],wind) ##selected vars

###CALIB AREAS
#Calibration Areas (50km buffers)
BuffAmer<-readOGR("Datos", layer = "BuffAM50km")
BuffMedit<-readOGR("Datos", layer = "BuffMED50km")
BuffSuraf<-readOGR("Datos", layer = "BuffSUD50km")


###BACKGROUNDS OF EACH STUDY AREAS DEFINED BY THE 50K BUFFERS
bg_inv<-mask(SelVars,BuffAmer)
bg_inv<-crop(bg_inv,BuffAmer)
plot(bg_inv[[1]])
points(amer[,c(3,2)])
inv_env <- na.omit(as.data.frame(bg_inv, xy=T))

bg_MED<-mask(SelVars,BuffMedit)
bg_MED<-crop(bg_MED,BuffMedit)
plot(bg_MED[[1]])
points(medit[,c(3,2)])
med_env <- na.omit(as.data.frame(bg_MED, xy=T))

bg_SUD<-mask(SelVars,BuffSuraf)
bg_SUD<-crop(bg_SUD,BuffSuraf)
plot(bg_SUD[[1]])
points(suraf[,c(3,2)])
sud_env <- na.omit(as.data.frame(bg_SUD, xy=T))

globalENV<-rbind(inv_env,med_env,sud_env)

####DATA FROM LOCALITIES
#Extracting data for the niche dynamics analysis in the invaded range
inv <- read.csv("Datos/Occs/todosamerica.csv")
inv<- cbind(inv[,c(3,2)], raster::extract(SelVars,inv[,c(3,2)]))
inv<-na.omit(inv)

#Extracting data for the niche dynamics analysis in the native ranges
nat_med <- read.csv("Datos/Occs/todosmediterraneo.csv")
nat_med<- cbind(nat_med[,c(3,2)], raster::extract(SelVars,nat_med[,c(3,2)]))
nat_med<-na.omit(nat_med)

nat_sud <- read.csv("Datos/Occs/todossurafrica.csv")
nat_sud<- cbind(nat_sud[,c(3,2)], raster::extract(SelVars,nat_sud[,c(3,2)]))
nat_sud<-na.omit(nat_sud)

###ENV-PCA
#The PCA is calibrated on all the sites of the study area, including both native and invaded ranges 
#(same as PCAenv in Broenniman et al. 2012)
pca.env <- dudi.pca(globalENV[,3:8],scannf=F,nf=2)

#Plot Variables Contribution with ecospat.plot.contrib()
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)

#Predict the scores on the axes
# PCA scores for the whole study area
scores.globclim <- pca.env$li

# PCA scores for the species native distribution
scores.sp.nat_med <- suprow(pca.env,nat_med[,3:8])$li
scores.sp.nat_sud <- suprow(pca.env,nat_sud[,3:8])$li

# PCA scores for the species invasive distribution
scores.sp.inv <- suprow(pca.env,inv[,3:8])$li

# PCA scores for the whole native study area
scores.clim.nat_med <- suprow(pca.env,med_env[,3:8])$li
scores.clim.nat_sud <- suprow(pca.env,sud_env[,3:8])$li

# PCA scores for the whole invaded study area
scores.clim.inv <- suprow(pca.env,inv_env[,3:8])$li

#Con convex hulls
#scores.clim.inv1 <- suprow(pca.env,inv_env1[,3:8])$li
#scores.clim.inv2 <- suprow(pca.env,inv_env2[,3:8])$li
#scores.clim.inv<-rbind(scores.clim.inv1,scores.clim.inv2)

# gridding the native niche
grid.clim.nat_med <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.nat_med,
                                       sp=scores.sp.nat_med, R=250,
                                       th.sp=0)

grid.clim.nat_sud <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.nat_sud,
                                       sp=scores.sp.nat_sud, R=250,
                                       th.sp=0)

# gridding the invasive niche MED
grid.clim.inv <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.inv,
                                       sp=scores.sp.inv, R=250,
                                       th.sp=0)

#2.5.3 Calculate Niche Overlap with ecospat.niche.overlap()
# Compute Schoener's D, index of niche overlap
D.overlap_med <- ecospat.niche.overlap (grid.clim.nat_med, grid.clim.inv, cor = TRUE)$D
D.overlap_med

D.overlap_sud <- ecospat.niche.overlap (grid.clim.nat_sud, grid.clim.inv, cor = TRUE)$D
D.overlap_sud

I.overlap_med <- ecospat.niche.overlap (grid.clim.nat_med, grid.clim.inv, cor = TRUE)$I
I.overlap_med

I.overlap_sud <- ecospat.niche.overlap (grid.clim.nat_sud, grid.clim.inv, cor = TRUE)$I
I.overlap_sud

# Perform the Niche Equivalency Test (according to Warren et al. (2008))
# It is reccomended to use at least 1000 replications for the equivalency test. 
eq.test_med <- ecospat.niche.equivalency.test(grid.clim.nat_med, grid.clim.inv,
                                          rep=1000, alternative = "greater")
ecospat.plot.overlap.test(eq.test_med, "D", "Equivalency MED-AMER")


eq.test_sud <- ecospat.niche.equivalency.test(grid.clim.nat_sud, grid.clim.inv,
                                              rep=1000, alternative = "greater")
ecospat.plot.overlap.test(eq.test_sud, "D", "Equivalency SUD-AMER")

# Niche Similarity Test
# Shifts randomly on niche (here the invasive niche) in the study area It is recomended to use at least 1000 replications
sim.test_med <- ecospat.niche.similarity.test(grid.clim.nat_med, grid.clim.inv,
                                          rep=1000, alternative = "greater",
                                          rand.type=2)
ecospat.plot.overlap.test(sim.test_med, "D", "Similarity MED-AMER")

sim.test_sud <- ecospat.niche.similarity.test(grid.clim.nat_sud, grid.clim.inv,
                                              rep=1000, alternative = "greater",
                                              rand.type=2)
ecospat.plot.overlap.test(sim.test_sud, "D", "Similarity SUD-AMER")

#Delimiting niche categories and quantifying niche dynamics in analogue climates
niche.dyn_med <- ecospat.niche.dyn.index (grid.clim.nat_med, grid.clim.inv, intersection = 0.1)
niche.dyn_sud <- ecospat.niche.dyn.index (grid.clim.nat_sud, grid.clim.inv, intersection = 0.1)


#Visualizing niche categories, niche dynamics and climate analogy between ranges
#with ecospat.plot.niche.dyn() Plot niche overlap

### MEDITERRANEO-AMERICA
par(mfrow=c(2,2))
ecospat.plot.niche.dyn(grid.clim.nat_med, grid.clim.inv, title= "Niche Overlap Mediterranean",quant=0.1,  
                       name.axis1="PC1",name.axis2="PC2",colz1 = c("royalblue3", alpha=0.15), colZ1="royalblue1", 
                       colZ2 = "red3",interest=2)
ecospat.shift.centroids(scores.sp.nat_med, scores.sp.inv, scores.clim.nat_med, scores.clim.inv)

### SOUTHAFRICA-AMERICA
ecospat.plot.niche.dyn(grid.clim.nat_sud, grid.clim.inv, title= "Niche Overlap SouthAfrica", quant=0.1,  
                       name.axis1="PC1",name.axis2="PC2",colz1 = c("royalblue3", alpha=0.15), colZ1="royalblue1", 
                       colZ2 = "red3",interest=2)
ecospat.shift.centroids(scores.sp.nat_sud, scores.sp.inv, scores.clim.nat_sud, scores.clim.inv)

