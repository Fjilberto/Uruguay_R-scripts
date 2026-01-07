#FJ 2025

rm(list=ls()) # borra el ambiente
graphics.off() # Limpiar la lista de graficos

library(tess3r)
library(LEA)
library(tidyverse)
library(ggspatial)
library(assignPOP)
library(vegan)
library(adegenet)

#Set files folder
setwd("D:/Trabajo/Paper/uruguay/ambiental/")
setwd("/media/oxus/LS/Trabajo/Paper/uruguay/ambiental")

#---- DATA MANIPULACION ----

### Convert structure to geno and lfml
struct2geno("uruprocesado80v2_allref_id_pop.str",ploidy = 2,FORMAT = 2,extra.row = 1,extra.column = 2)

### Algunas funciones exigen no tener NA en el archivo de genotipos y piden imputarlos

# pca
pc = pca("uruprocesado80v2_allref_id_pop.str.lfmm", scale = TRUE)
tw = tracy.widom(pc)

# plot the percentage of variance explained by each component
plot(tw$percentage, pch = 19, col = "darkblue", cex = .8)

project = NULL
project = snmf("uruprocesado80v2_allref_id_pop.str.geno", 
               K = 1:10, entropy = TRUE, repetitions = 10, project = "new")

# plot cross-entropy criterion for all runs in the snmf project
plot(project, col = "blue", pch = 19, cex = 1.2)

# select the best run for K = 4 clusters

best = which.min(cross.entropy(project, K = 3))
my.colors <- c("tomato", "lightblue",
               "olivedrab", "gold","dodgerblue","plum")
barchart(project, K = 3, run = best,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),labels = bp$order, las=1,cex.axis = .4)

project.missing = snmf("uruprocesado80v2_allref_id_pop.str.lfmm", K = 3,
                       entropy = TRUE, repetitions = 10,
                       project = "new")

impute(project.missing, "uruprocesado80v2_allref_id_pop.str.lfmm",
       method = 'mode', K = 3, run = best)

#### clonar coordenadas y datos ambientales

# coordenadas

structure<-read.csv("uruprocesado80v2_allref_id_pop.str",sep = "\t")
conteo<-as.data.frame(table(structure$pop))
conteo$ind<-conteo$Freq/2;sum(conteo$ind)

coordenadas<-read.csv("pointfile2.txt",sep = ",",dec = ".")
coordenadas$loc <- as.factor(coordenadas$loc)
coordenadas$loc <- factor(coordenadas$loc,
                          levels=c("PBUY","LPUY","PEUY","PRUY","BAAR","MPAR","BBAR","FKUK","BACL","PMAR",
                                   "CTCL","RMCL","PCCL","IPCL",
                                   "SGCL","BLCL","PACL","ICCL","PPCL","BICL","CMCL",
                                   "CAAR","ALAR","PNCL","PWCL"))
coordenadas <- coordenadas[order(coordenadas$loc), ]
rownames(coordenadas) <- NULL
coordenadas<-coordenadas[-18,]
coordenadas<-coordenadas[,-1]

if (nrow(coordenadas) != length(conteo$ind)) {
  stop("El número de filas del dataframe debe ser igual a la longitud del vector de copias.")
}

# Crear una lista para almacenar las filas duplicadas
filas_duplicadas <- list()

# Iterar a través de cada fila del dataframe y el número correspondiente del vector
for (i in 1:nrow(coordenadas)) {
  fila_a_duplicar <- coordenadas[i, ]  # Seleccionar la fila actual
  num_repeticiones <- conteo$ind[i] # Obtener el número de veces a copiar
  
  # Duplicar la fila y agregarla a la lista
  for (j in 1:num_repeticiones) {
    filas_duplicadas[[length(filas_duplicadas) + 1]] <- fila_a_duplicar
  }
}

# Convertir la lista de filas duplicadas a un nuevo dataframe
nuevo_coordenadas <- do.call(rbind, filas_duplicadas)
write.table(nuevo_coordenadas,file = "nuevo_coordenadas.txt",sep = ",",dec = ".",
            row.names = F,col.names = F,quote = F)

nuevo_coordenadas2<-as.matrix(nuevo_coordenadas)
nuevo_coordenadas2[, c(1, 2)] <- nuevo_coordenadas2[, c(2, 1)]

plot(nuevo_coordenadas2, pch = 19, cex = .5, 
     xlab = "Longitude (°E)", ylab = "Latitude (°N)")
maps::map(add = T, interior = T)

## datos ambientales

ambientales<-read.csv("data_ambiental.env",sep = ",",dec = ".",header = T)

# Verificar que el número de filas del dataframe coincida con la longitud del vector
if (nrow(ambientales) != length(conteo$ind)) {
  stop("El número de filas del dataframe debe ser igual a la longitud del vector de copias.")
}

# Crear una lista para almacenar las filas duplicadas
filas_duplicadas <- list()

# Iterar a través de cada fila del dataframe y el número correspondiente del vector
for (i in 1:nrow(ambientales)) {
  fila_a_duplicar <- ambientales[i, ]  # Seleccionar la fila actual
  num_repeticiones <- conteo$ind[i] # Obtener el número de veces a copiar
  
  # Duplicar la fila y agregarla a la lista
  for (j in 1:num_repeticiones) {
    filas_duplicadas[[length(filas_duplicadas) + 1]] <- fila_a_duplicar
  }
}

# Convertir la lista de filas duplicadas a un nuevo dataframe
nuevo_ambientales <- do.call(rbind, filas_duplicadas)
write.table(nuevo_ambientales,file = "nuevo_ambientales.env",sep = ",",dec = ".",
            row.names = F,col.names = F,quote = F)

#---- Tess3: Analisis geoespacial ----

X<-read.csv("uruprocesado80v2_allref_id_pop.str.lfmm_imputed.lfmm",sep = " ",header = F)

tess3.obj <- tess3(X = X, coord = nuevo_coordenadas2, K = 1:8, 
                   method = "projected.ls", ploidy = 2, openMP.core.num = 4) 

plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

map_bathy <- marmap::getNOAA.bathy(lon1 = -80, lon2 = -50, lat1 = -30,
                                   lat2 = -60, res=10,
                                   keep=T)

map_bathy1 = - map_bathy

## depth control: exclude areas that are deeper than 4000 meters 
## (in the end, we did not use this feature)
for (i in 1:180) {
  for (j in 1:180) {
    if (map_bathy1[i,j] > 100 | map_bathy1[i,j] < 0) {
      map_bathy1[i,j] <- -1  
    }
  }
}

summary(map_bathy1)

## convert the bathymetry map to a raster

asc_raster <- marmap::as.raster(map_bathy1)

## Write raster to file for later use

raster::writeRaster(asc_raster, "D:/Trabajo/Paper/uruguay/ambiental/uru.asc", 
                    overwrite=TRUE)

# retrieve tess3 Q matrix for K cluster 
q.matrix2 <- qmatrix(tess3.obj, K = 2)
q.matrix3 <- qmatrix(tess3.obj, K = 3)
q.matrix4 <- qmatrix(tess3.obj, K = 4)
q.matrix7 <- qmatrix(tess3.obj, K = 7)

# STRUCTURE-like barplot for the Q-matrix 
barplot(q.matrix2, border = NA, space = 0, 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "Ancestry matrix") -> bp
## Use CreatePalette() to define color palettes.
axis(1, at = 1:nrow(q.matrix2), labels = bp$order, las = 3, cex.axis = .4)

## Test plot of the perm-dataset with adapted raster and K = 4

plot(q.matrix2, nuevo_coordenadas2, method="map.max", cex=1,
     raster.filename = "D:/Trabajo/Paper/uruguay/ambiental/uru.asc", 
     interpol=FieldsKrigModel(9.5),
     main=paste0("Ancestry coefficients"), resolution=c(300,300),
     xlab="longitude", ylab="Latitude")

## I assembled the final plot with ggplot, and the ggtess3Q function from the tess3R package.
## first, I called a costum color palette. Also, I need a world map as a backdrop for the plot,
## which can be created using the rworldmap package. 

## Get a map
testmap <- rworldmap::getMap(resolution = "high", projection=NA)

## create a color palette and shades of the color palette with CreatePalette
my.colors <- c("#1D72F5","#DF0101","#77CE61","#FF9326","#A945FF",
               "#0089B2","#FDF060","#FFA6B2")
my.palette <- CreatePalette(my.colors, 15)

mycol=c("#E9CE53","#EC7B74") #k2
mycol=c("#56B9D7","#EC7B74","#E9CE53") #k3
mycol=c("#56B9D7","#E9CE53","#D692C6","#EC7B74") #k4
mycol=c("#B0B0B0","#EC7B74","#E9CE53","#56B9D7","#D692C6","#79B195","#B19470") #K7

## test plot for K = 4, perm-dataset
pl <- ggtess3Q(q.matrix7, nuevo_coordenadas2, 
               raster.filename = "D:/Trabajo/Paper/uruguay/ambiental/uru.asc",
               col.palette = mycol)
pl4 <- pl + geom_path(data = testmap, aes(x = long, y = lat, group = group)) +
  xlim(-78, -52) + 
  ylim(-58, -33) + 
  coord_equal()  +
  geom_point(data = as.data.frame(nuevo_coordenadas2), aes(x = new_lat, y = new_lon), size = 3, 
             color = "black") + 
  xlab("Longitute") +
  ylab("Latitude") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_colour_continuous(mycol) #+
  #annotation_scale(location = "bl", width_hint = 0.15) 

pl4

ggsave("tess3_k7.tiff", plot = last_plot(),
       scale = 2, width = 21.94, height = 10.59 , 
       units = "cm", dpi = 300)

#---- LEA Population differentiation ----

# Genome scan for selection: population differentiation tests
p = snmf.pvalues(project,
                 entropy = TRUE,
                 ploidy = 2,
                 K = 3)
pvalues = p$pvalues
par(mfrow = c(2,1))
hist(pvalues, col = "orange")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .5)

######### Data ambiental

# Ecological association tests using lfmm

Y<-read.csv("uruprocesado80v2_allref_id_pop.str.lfmm_imputed.lfmm",sep = " ")
X<-read.csv("nuevo_ambientales.env",sep = ",",dec = ".")

mod.lfmm2 <- lfmm2(input = Y, env = X[,-12], K = 4)
plot(mod.lfmm2@U, col = "grey", pch = 19,
     xlab = "Factor 1",
     ylab = "Factor 2")

pv <- lfmm2.test(object = mod.lfmm2,
                 input = Y,
                 env = X[,-12],
                 full = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.1/510), lty = 2, col = "orange")

#---- AssignPOP ----

# genetic data

popnames<-c("LPUY","PEUY","PRUY","BSAR","BBAR","MPAR","FKUK","BACL","PMAR",
            "CTCL","RMCL","PCCL","IPCL",
            "SGCL","BLCL","PACL","PPCL","BICL","CMCL",
            "CAAR","ALAR","PNCL","PWCL")

gen_data<-read.Genepop( "uruprocesado80v2_onlystudy_nopbuy.gen", pop.names=popnames, haploid = FALSE)

# ambiental data

amb_data<-read.csv("nuevo_ambientales.env",sep = ",",dec = ".",header = F)
amb_data<-amb_data[-c(1:4),]
datos_amb <- read.csv("data_ambiental.txt", row.names = 1)
colnames(amb_data)<-c(names(datos_amb))

pop_label<-stringr::str_sub(gen_data$SampleID,start = 1, end = 4)
pop_label[pop_label=="ICCL"]<-"PPCL"
ambpop_data <- cbind(gen_data$SampleID, amb_data, pop_label)

names(ambpop_data)[1]<-"ID"
ambpop_data$pop_label<-as.factor(ambpop_data$pop_label)

# diferentes

n_filas <- nrow(ambpop_data)
secuencia_a_sumar <- seq(0.00011, by = 0.00011, length.out = n_filas)

# Sumar la secuencia a todas las columnas numéricas
ambpop_data <- ambpop_data %>%
  mutate(across(where(is.numeric), function(x) x + secuencia_a_sumar))

# genetic + ambiental data

id_amb<-cbind(gen_data$SampleID,amb_data)
write.csv(id_amb,file = "id_amb.csv",row.names = F)

gen_amb_data<-compile.data(gen_data,"id_amb.csv")

#### Resampling cross-validation

#Monte_carlo

#the code below performs Monte-Carlo cross-validation, with using 50%, 70%, and 90% of 
#random individuals from each population (arg. train.inds) crossed by top 10%, 25%, 50% of high FST loci,
#and all loci (arg. train.loci, loci.sample="fst") as training data. Each combination of 
#training data is resampled 30 times (arg. iterations). 

#As a result, it performs a total of 360 assignment tests 
#(3 levels of training individuals by 4 levels of training loci by 30 iterations).

#only genetic
assign.MC( gen_data, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="svm", dir="Result_MC_genetic/")

#only ambiental
assign.MC( ambpop_data, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="svm", dir="Result_MC_amb/",scaled = F)

#genetic + ambiental
assign.MC( gen_amb_data, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="svm", dir="Result_MC_gen_amb/",scaled = F)

#k-fold cross validation

#the code divides individuals from each population into k groups (or folds). 
#In each fold, 10%, 25%, and 50% of random loci (loci.sample = "random") 
#as well as all loci were sampled as training data. As a result, 
#the K-fold cross-validation performs a total of x assignment tests 
#(for 3 groups= 48 = (3+4+5) folds * 4 levels of training loci).

#only genetic
assign.kfold( gen_data, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="lda", dir="Result_kfold_genetic/")

#only ambiental
assign.kfold( ambpop_data, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="lda", dir="Result_kfold_ambiental/")

#genetic + ambiental
assign.kfold( gen_amb_data, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="lda", dir="Result_kfold_gen_amb/")

## Assgnment accuracy

accu_kf_gen<-accuracy.kfold(dir="Result_kfold_genetic/")
accuracy.plot(accu_kf_gen, pop = "all")

accuracy.plot(accu_kf_gen, pop=c("LPUY","PEUY","PRUY","BSAR","BBAR","MPAR")) +
  ylim(0, 1) + #Set y limit between 0 and 1
  annotate("segment",x=0.4,xend=3.6,y=0.33,yend=0.33,colour="red",size=1) + #Add a red horizontal line at y = 0.33 (null assignment rate for 3 populations)
  ggtitle("Monte-Carlo cross-validation using genetic loci")+ #Add a plot title
  theme(plot.title = element_text(size=20, face="bold")) #Edit plot title text size

## informative loci

check.loci(dir = "Result_MC_genetic/", top.loci = 20)

##membership plot only for kfold

membership.plot(dir = "Result_kfold_genetic/")

## Assigment matrix

assign.matrix( dir="Result_kfold_genetic/")

####### k2

popnames<-c("chilensis_like","platensis_like")

gen_data2<-read.Genepop( "uruprocesado80v2_onlystudy_nopbuy_k2.gen", pop.names=popnames, haploid = FALSE)

pop_label<-gen_data2$DataMatrix$popNames_vector
ambpop_data <- cbind(gen_data2$SampleID, amb_data, pop_label)

names(ambpop_data)[1]<-"ID"
ambpop_data$pop_label<-as.factor(ambpop_data$pop_label)
str(ambpop_data)

# diferentes

n_filas <- nrow(ambpop_data)
secuencia_a_sumar <- seq(0.00011, by = 0.00011, length.out = n_filas)

# Sumar la secuencia a todas las columnas numéricas
ambpop_data <- ambpop_data %>%
  mutate(across(where(is.numeric), function(x) x + secuencia_a_sumar))

# genetic + ambiental data

id_amb<-cbind(gen_data2$SampleID,amb_data)
write.csv(id_amb,file = "id_amb2.csv",row.names = F)

gen_amb_data<-compile.data(gen_data2,"id_amb2.csv")

## Monte Carlo

assign.MC( gen_data2, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="lda", dir="Result_MC_genetic2_fst/")

assign.MC( gen_data2, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="svm", dir="Result_MC_genetic2_fst_svm/")

assign.MC( gen_data2, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="tree", dir="Result_MC_genetic2_fst_tree/")

assign.MC( gen_amb_data, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="lda", dir="Result_MC_genamb2_fst/")

assign.MC( gen_amb_data, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="svm", dir="Result_MC_genamb2_fst_svm/")

assign.MC( gen_amb_data, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="tree", dir="Result_MC_genamb2_fst_tree/")

accu_mc_gen2<-accuracy.MC(dir="Result_MC_genetic2_fst/")
accu_mc_gen2<-accuracy.MC(dir="Result_MC_genetic2_fst_svm/")
accu_mc_gen2<-accuracy.MC(dir="Result_MC_genetic2_fst_tree/")

accu_mc_gen2<-accuracy.MC(dir="Result_MC_genamb2_fst/")
accu_mc_gen2<-accuracy.MC(dir="Result_MC_genamb2_fst_svm/")
accu_mc_gen2<-accuracy.MC(dir="Result_MC_genamb2_fst_tree/")

accuracy.plot(accu_mc_gen2, pop = "all")

accuracy.plot(accu_mc_gen2, pop=c("chilensis_like","platensis_like")) +
  ylim(0, 1) + #Set y limit between 0 and 1
  annotate("segment",x=0.4,xend=3.6,y=0.33,yend=0.33,colour="red",size=1) + #Add a red horizontal line at y = 0.33 (null assignment rate for 3 populations)
  ggtitle("Monte-Carlo cross-validation using genetic loci")+ #Add a plot title
  theme(plot.title = element_text(size=20, face="bold")) 

## Assigment matrix

assign.matrix( dir="Result_MC_genetic2_fst/")
assign.matrix( dir="Result_MC_genetic2_fst_svm/")
assign.matrix( dir="Result_MC_genetic2_fst_tree/")

assign.matrix( dir="Result_MC_genamb2_fst/")
assign.matrix( dir="Result_MC_genamb2_fst_svm/")
assign.matrix( dir="Result_MC_genamb2_fst_tree/")

## Kfold

assign.kfold( gen_data2, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="lda", dir="Result_kfold_genetic2_fst/")

assign.kfold( gen_data2, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="svm", dir="Result_kfold_genetic2_fst_svm/")

assign.kfold( gen_data2, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="tree", dir="Result_kfold_genetic2_fst_tree/")

assign.kfold( gen_amb_data, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="lda", dir="Result_kfold_genamb2_fst/")

assign.kfold( gen_amb_data, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="svm", dir="Result_kfold_genamb2_fst_svm/")

assign.kfold( gen_amb_data, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="tree", dir="Result_kfold_genamb2_fst_tree/")


accu_kf_gen2<-accuracy.kfold(dir="Result_kfold_genetic2_fst/")
accu_kf_gen2<-accuracy.kfold(dir="Result_kfold_genetic2_fst_svm/")
accu_kf_gen2<-accuracy.kfold(dir="Result_kfold_genetic2_fst_tree/")

accu_kf_gen2<-accuracy.kfold(dir="Result_kfold_genamb2_fst/")
accu_kf_gen2<-accuracy.kfold(dir="Result_kfold_genamb2_fst_svm/")
accu_kf_gen2<-accuracy.kfold(dir="Result_kfold_genamb2_fst_tree/")

accuracy.plot(accu_kf_gen2, pop = "all")

accuracy.plot(accu_kf_gen2, pop=c("chilensis_like","platensis_like")) +
  ylim(0, 1) + #Set y limit between 0 and 1
  annotate("segment",x=0.4,xend=3.6,y=0.33,yend=0.33,colour="red",size=1) + #Add a red horizontal line at y = 0.33 (null assignment rate for 3 populations)
  ggtitle("Monte-Carlo cross-validation using genetic loci")+ #Add a plot title
  theme(plot.title = element_text(size=20, face="bold")) 

##membership plot only for kfold

membership.plot(dir = "Result_kfold_genetic2_fst/")
membership.plot(dir = "Result_kfold_genetic2_fst_svm/")
membership.plot(dir = "Result_kfold_genetic2_fst_tree/")

membership.plot(dir = "Result_kfold_genamb2_fst/")
membership.plot(dir = "Result_kfold_genamb2_fst_svm/")
membership.plot(dir = "Result_kfold_genamb2_fst_tree/")

## Assigment matrix

assign.matrix( dir="Result_kfold_genetic2_fst/",k.fold = 2)
assign.matrix( dir="Result_kfold_genetic2_fst_svm/",k.fold = 2)
assign.matrix( dir="Result_kfold_genetic2_fst_tree/",k.fold = 2)

assign.matrix( dir="Result_kfold_genamb2_fst/",k.fold = 2)
assign.matrix( dir="Result_kfold_genamb2_fst_svm/",k.fold = 2)
assign.matrix( dir="Result_kfold_genamb2_fst_tree/",k.fold = 2)

####### k3
popnames<-c("uruguay","chilensis_like","platensis_like")

gen_data3<-read.Genepop( "uruprocesado80v2_onlystudy_nopbuy_k3.gen", pop.names=popnames, haploid = FALSE)

pop_label<-gen_data3$DataMatrix$popNames_vector
ambpop_data <- cbind(gen_data3$SampleID, amb_data, pop_label)

names(ambpop_data)[1]<-"ID"
ambpop_data$pop_label<-as.factor(ambpop_data$pop_label)
str(ambpop_data)

# diferentes

n_filas <- nrow(ambpop_data)
secuencia_a_sumar <- seq(0.00011, by = 0.00011, length.out = n_filas)

# Sumar la secuencia a todas las columnas numéricas
ambpop_data <- ambpop_data %>%
  mutate(across(where(is.numeric), function(x) x + secuencia_a_sumar))

# genetic + ambiental data

id_amb<-cbind(gen_data3$SampleID,amb_data)
write.csv(id_amb,file = "id_amb3.csv",row.names = F)

gen_amb_data<-compile.data(gen_data3,"id_amb3.csv")

#Montecarlo

assign.MC( gen_data3, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="lda", dir="Result_MC_genetic3_fst/")

assign.MC( gen_data3, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="svm", dir="Result_MC_genetic3_fst_svm/")

assign.MC( gen_data3, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="tree", dir="Result_MC_genetic3_fst_tree/")

assign.MC( gen_amb_data, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="lda", dir="Result_MC_genamb3_fst/")

assign.MC( gen_amb_data, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="svm", dir="Result_MC_genamb3_fst_svm/")

assign.MC( gen_amb_data, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="tree", dir="Result_MC_genamb3_fst_tree/")

accu_mc_gen3<-accuracy.MC(dir="Result_MC_genetic3_fst/")
accu_mc_gen3<-accuracy.MC(dir="Result_MC_genetic3_fst_svm/")
accu_mc_gen3<-accuracy.MC(dir="Result_MC_genetic3_fst_tree/")

accuracy.plot(accu_mc_gen3, pop = "all")

accuracy.plot(accu_mc_gen3, pop=c("uruguay","chilensis_like","platensis_like")) +
  ylim(0, 1) + #Set y limit between 0 and 1
  annotate("segment",x=0.4,xend=3.6,y=0.33,yend=0.33,colour="red",size=1) + #Add a red horizontal line at y = 0.33 (null assignment rate for 3 populations)
  ggtitle("Monte-Carlo cross-validation using genetic loci")+ #Add a plot title
  theme(plot.title = element_text(size=20, face="bold")) 

## Assigment matrix

assign.matrix( dir="Result_MC_genetic3_fst/")
assign.matrix( dir="Result_MC_genetic3_fst_svm/")
assign.matrix( dir="Result_MC_genetic3_fst_tree/")

assign.matrix( dir="Result_MC_genamb3_fst/")
assign.matrix( dir="Result_MC_genamb3_fst_svm/")
assign.matrix( dir="Result_MC_genamb3_fst_tree/")

#kfold

assign.kfold( gen_data3, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="lda", dir="Result_kfold_genetic3_fst/")

assign.kfold( gen_data3, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="svm", dir="Result_kfold_genetic3_fst_svm/")

assign.kfold( gen_data3, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="tree", dir="Result_kfold_genetic3_fst_tree/")

assign.kfold( gen_amb_data, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="lda", dir="Result_kfold_genamb3_fst/")

assign.kfold( gen_amb_data, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="svm", dir="Result_kfold_genamb3_fst_svm/")

assign.kfold( gen_amb_data, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="tree", dir="Result_kfold_genamb3_fst_tree/")

accu_kf_gen3<-accuracy.kfold(dir="Result_kfold_genetic3_fst/")
accu_kf_gen3<-accuracy.kfold(dir="Result_kfold_genetic3_fst_svm/")
accu_kf_gen3<-accuracy.kfold(dir="Result_kfold_genetic3_fst_tree/")
accuracy.plot(accu_kf_gen3, pop = "all")


accuracy.plot(accu_kf_gen3, pop=c("uruguay","chilensis_like","platensis_like")) +
  ylim(0, 1) + #Set y limit between 0 and 1
  annotate("segment",x=0.4,xend=3.6,y=0.33,yend=0.33,colour="red",size=1) + #Add a red horizontal line at y = 0.33 (null assignment rate for 3 populations)
  ggtitle("Monte-Carlo cross-validation using genetic loci")+ #Add a plot title
  theme(plot.title = element_text(size=20, face="bold")) 

##membership plot only for kfold

membership.plot(dir = "Result_kfold_genetic2_fst/")
membership.plot(dir = "Result_kfold_genetic2_fst_svm/")
membership.plot(dir = "Result_kfold_genetic2_fst_tree/")

membership.plot(dir = "Result_kfold_genamb3_fst/")
membership.plot(dir = "Result_kfold_genamb3_fst_svm/")
membership.plot(dir = "Result_kfold_genamb3_fst_tree/")

## Assigment matrix

assign.matrix( dir="Result_kfold_genetic3_fst/",k.fold = 3)
assign.matrix( dir="Result_kfold_genetic3_fst_svm/",k.fold = 3)
assign.matrix( dir="Result_kfold_genetic3_fst_tree/",k.fold = 3)

assign.matrix( dir="Result_kfold_genamb3_fst/",k.fold = 3)
assign.matrix( dir="Result_kfold_genamb3_fst_svm/",k.fold = 3)
assign.matrix( dir="Result_kfold_genamb3_fst_tree/",k.fold = 3)

####### k4
popnames<-c("uruguay","chilensis_like","hybrid_zone","platensis_like")

gen_data4<-read.Genepop( "uruprocesado80v2_onlystudy_nopbuy_k4.gen", pop.names=popnames, haploid = FALSE)

pop_label<-gen_data4$DataMatrix$popNames_vector
ambpop_data <- cbind(gen_data4$SampleID, amb_data, pop_label)

names(ambpop_data)[1]<-"ID"
ambpop_data$pop_label<-as.factor(ambpop_data$pop_label)
str(ambpop_data)

# diferentes

n_filas <- nrow(ambpop_data)
secuencia_a_sumar <- seq(0.00011, by = 0.00011, length.out = n_filas)

# Sumar la secuencia a todas las columnas numéricas
ambpop_data <- ambpop_data %>%
  mutate(across(where(is.numeric), function(x) x + secuencia_a_sumar))

# genetic + ambiental data

id_amb<-cbind(gen_data4$SampleID,amb_data)
write.csv(id_amb,file = "id_amb4.csv",row.names = F)

gen_amb_data<-compile.data(gen_data4,"id_amb4.csv")

#Montecarlo

assign.MC( gen_data4, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="lda", dir="Result_MC_genetic4_fst/")

assign.MC( gen_data4, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="svm", dir="Result_MC_genetic4_fst_svm/")

assign.MC( gen_data4, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="tree", dir="Result_MC_genetic4_fst_tree/")

assign.MC( gen_amb_data, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="lda", dir="Result_MC_genamb4_fst/")

assign.MC( gen_amb_data, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="svm", dir="Result_MC_genamb4_fst_svm/")

assign.MC( gen_amb_data, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="tree", dir="Result_MC_genamb4_fst_tree/")

accu_mc_gen4<-accuracy.MC(dir="Result_MC_genetic4_fst/")
accu_mc_gen4<-accuracy.MC(dir="Result_MC_genetic4_fst_svm/")
accu_mc_gen4<-accuracy.MC(dir="Result_MC_genetic4_fst_tree/")

accuracy.plot(accu_mc_gen4, pop = "all")

accuracy.plot(accu_mc_gen4, pop=c("uruguay","chilensis_like","platensis_like")) +
  ylim(0, 1) + #Set y limit between 0 and 1
  annotate("segment",x=0.4,xend=3.6,y=0.33,yend=0.33,colour="red",size=1) + #Add a red horizontal line at y = 0.33 (null assignment rate for 3 populations)
  ggtitle("Monte-Carlo cross-validation using genetic loci")+ #Add a plot title
  theme(plot.title = element_text(size=20, face="bold")) 

## Assigment matrix

assign.matrix( dir="Result_MC_genetic4_fst/")
assign.matrix( dir="Result_MC_genetic4_fst_svm/")
assign.matrix( dir="Result_MC_genetic4_fst_tree/")

assign.matrix( dir="Result_MC_genamb4_fst/")
assign.matrix( dir="Result_MC_genamb4_fst_svm/")
assign.matrix( dir="Result_MC_genamb4_fst_tree/")

##kfold

assign.kfold( gen_data4, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="lda", dir="Result_kfold_genetic4_fst/")

assign.kfold( gen_data4, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="svm", dir="Result_kfold_genetic4_fst_svm/")

assign.kfold( gen_data4, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="tree", dir="Result_kfold_genetic4_fst_tree/")

assign.kfold( gen_amb_data, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="lda", dir="Result_kfold_genamb4_fst/")

assign.kfold( gen_amb_data, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="svm", dir="Result_kfold_genamb4_fst_svm/")

assign.kfold( gen_amb_data, k.fold=c(2,3,4,5,6,7,8), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="fst", model="tree", dir="Result_kfold_genamb4_fst_tree/")

accu_kf_gen4<-accuracy.kfold(dir="Result_kfold_genetic4_fst/")
accu_kf_gen4<-accuracy.kfold(dir="Result_kfold_genetic4_fst_svm/")
accu_kf_gen4<-accuracy.kfold(dir="Result_kfold_genetic4_fst_tree/")
accuracy.plot(accu_kf_gen4, pop = "all")

accuracy.plot(accu_kf_gen4, pop=c("uruguay","chilensis_like","hybrid_zone","platensis_like")) +
  ylim(0, 1) + #Set y limit between 0 and 1
  annotate("segment",x=0.4,xend=3.6,y=0.33,yend=0.33,colour="red",size=1) + #Add a red horizontal line at y = 0.33 (null assignment rate for 3 populations)
  ggtitle("Monte-Carlo cross-validation using genetic loci")+ #Add a plot title
  theme(plot.title = element_text(size=20, face="bold")) 

##membership plot only for kfold

membership.plot(dir = "Result_kfold_genetic4_fst/")
membership.plot(dir = "Result_kfold_genetic4_fst_svm/")
membership.plot(dir = "Result_kfold_genetic4_fst_tree/")

membership.plot(dir = "Result_kfold_genamb4_fst/")
membership.plot(dir = "Result_kfold_genamb4_fst_svm/")
membership.plot(dir = "Result_kfold_genamb4_fst_tree/")

## Assigment matrix

assign.matrix( dir="Result_kfold_genetic4_fst/")
assign.matrix( dir="Result_kfold_genetic4_fst_svm/")
assign.matrix( dir="Result_kfold_genetic4_fst_tree/")

assign.matrix( dir="Result_kfold_genamb4_fst/")
assign.matrix( dir="Result_kfold_genamb4_fst_svm/")
assign.matrix( dir="Result_kfold_genamb4_fst_tree/")

#---- Vegan RDA ----

# Cargar datos genéticos (ejemplo con formato genepop)
datos_gen <- read.genepop("uruprocesado80v2_onlystudy_nopbuy.gen", ncode = 3)

##### Identificar loci monomorficos
isPoly(datos_gen) %>% summary

#Eliminar monomorficos

poly_loci = names(which(isPoly(datos_gen) == TRUE))
datos_gen = datos_gen[loc = poly_loci]
isPoly(datos_gen) %>% summary

# Convertir a frecuencias alélicas
frec_alelicas <- genind2genpop(datos_gen) %>% tab()

# Cargar datos ambientales (debe tener las mismas muestras que los datos genéticos)
datos_amb <- read.csv("data_ambiental.txt", row.names = 1)
nombres_nuevos <- sub("^[^.]+\\.", "", names(datos_amb))
names(datos_amb) <- nombres_nuevos
datos_amb <- datos_amb[-c(3,19),]

rownames(frec_alelicas)<-rownames(datos_amb)

# Asegurarse de que las filas coincidan
frec_alelicas <- frec_alelicas[rownames(datos_amb), ]

# Estandarizar los datos ambientales (importante para RDA)
datos_amb_std <- scale(datos_amb)
datos_amb_std<-as.data.frame(datos_amb_std)

# Eliminar columnas con NA

datos_amb_std2 <- datos_amb_std %>%
  select(where(~ !anyNA(.)))

# Correlacion entre las variables ambientales

PerformanceAnalytics::chart.Correlation(datos_amb_std2)

library(ggcorrplot)
library(reshape2)

corr = round(cor(datos_amb_std2), 2)

# Get p-value matrix
p.df = as.data.frame(ggcorrplot::cor_pmat(datos_amb_std2))

# Function to get asteriks
labs.function = function(x){
  case_when(x >= 0.05 ~ "",
            x < 0.05 & x >= 0.01 ~ "*",
            x < 0.01 & x >= 0.001 ~ "**",
            x < 0.001 ~ "***")
}

# Get asteriks matrix based on p-values
p.labs = p.df  %>%                      
  mutate_all(labs.function)

# Reshaping asteriks matrix to match ggcorrplot data output
p.labs$Var1 = as.factor(rownames(p.labs))
p.labs = melt(p.labs, id.vars = "Var1", variable.name = "Var2", value.name = "lab")

# Initial ggcorrplot
cor.plot = ggcorrplot(corr, hc.order = TRUE, type = "lower",
                      lab = TRUE)

# Subsetting asteriks matrix to only those rows within ggcorrplot data
p.labs$in.df = ifelse(is.na(match(paste0(p.labs$Var1, p.labs$Var2), 
                                  paste0(cor.plot[["data"]]$Var1, cor.plot[["data"]]$Var2))),
                      "No", "Yes")

p.labs = select(filter(p.labs, in.df == "Yes"), -in.df)

# Add asteriks to ggcorrplot
cor.plot.labs = cor.plot + 
  geom_text(aes(x = p.labs$Var1, 
                y = p.labs$Var2), 
            label = p.labs$lab, 
            nudge_y = 0.25, 
            size = 5)
cor.plot.labs

ggsave("corr_plot.tiff", plot = last_plot(), dpi="print")

# Calcular matriz de correlación completa (SNPs × variables)
cor_matrix <- cor(frec_alelicas, datos_amb_std2,method = "spearman")

library(Hmisc)

datos_combinados <- cbind(frec_alelicas, datos_amb_std2)
resultados_cor <- rcorr(as.matrix(datos_combinados), type = "spearman")
p_valores <- resultados_cor$P
nombres_marcadores <- colnames(frec_alelicas)
nombres_ambientales <- colnames(datos_amb_std2)

p_valores_filtrados <- p_valores[nombres_marcadores, nombres_ambientales]
print(p_valores_filtrados)

# 1. Obtener la matriz de p-valores filtrada (como hicimos antes)
# Asumimos que `frec_alelicas` y `datos_amb_std2` son tus dataframes originales
library(Hmisc)
library(dplyr)
library(tidyr)
library(tibble)

datos_combinados <- cbind(frec_alelicas, datos_amb_std2)
resultados_cor <- rcorr(as.matrix(datos_combinados), type = "spearman")
p_valores_filtrados <- resultados_cor$P[colnames(frec_alelicas), colnames(datos_amb_std2)]

# 2. Transformar la matriz de p-valores a formato largo
p_valores_largos <- p_valores_filtrados %>%
  as.data.frame() %>%
  rownames_to_column("SNP") %>%
  pivot_longer(
    cols = -SNP,
    names_to = "Predictor",
    values_to = "p_value"
  )

# 3. Reorganizar la matriz de correlación al formato largo y unirla con los p-valores
resultados_finales <- resultados_cor$r[colnames(frec_alelicas), colnames(datos_amb_std2)] %>%
  as.data.frame() %>%
  rownames_to_column("SNP") %>%
  pivot_longer(
    cols = -SNP,
    names_to = "Predictor",
    values_to = "Correlation"
  ) %>%
  # Unir con la tabla de p-valores usando SNP y Predictor como claves
  left_join(p_valores_largos, by = c("SNP", "Predictor")) %>%
  mutate(
    N = 1:n()
  ) %>%
  # Eliminar la columna Eje_RDA
  select(N, SNP, Predictor, Correlation, p_value) %>%
  arrange(Predictor, desc(abs(Correlation))) %>%
  mutate(
    Correlation = round(Correlation, 2),
    p_value = round(p_value, 4) # Redondea los p-valores
  )

# Ver resultados
print(resultados_finales)

# Filtrar correlaciones con |r| > 0.8 (ajusta este umbral)
cor_extremas <- resultados_finales %>%
  filter(abs(p_value) < 0.05 & abs(Correlation) >= 0.8) %>%
  arrange(desc(abs(p_value)))

# Ver resultados
print(cor_extremas, n = Inf)  # Muestra todas las filas
write.table(cor_extremas,file = "cor_extremas.txt",sep = ";",dec = ",",row.names = F)

datos_gen_pop <- genind2genpop(datos_gen)
frecuencies_cor_candidates <- allele_frequencies[, unique(cor_extremas$SNP)]
write.table(frecuencies_cor_candidates,file = "frecuencias_cor.txt",sep = ";",dec = ",",row.names = T)

##### RDA

rda_result <- rda(frec_alelicas ~ OceanTemperature+Salinity+Nitrate+
                                  Phosphate+Silicate+DissolvedMolecularOxygen+
                                  pH+DissolvedIron+Chlorophyll+
                                  PhotosyntheticallyAvailableRadiation+TotalCloudFraction+AirTemperature+
                                  Precipitation+TidalHeight+TidalRange+
                                  EastwardSeaWaterVelocity+NorthwardSeaWaterVelocity, 
                  data = datos_amb_std2)

anova.cca(rda_result, permutations = 100000)
anova.cca(rda_result, permutations = 100000, by="terms")

# 2. Identificar SNPs significativos (p < 0.05)
signif_snps <- function(rda_obj, alpha = 0.05){
  # Calcular puntuaciones de los SNPs
  loadings <- scores(rda_obj, choices = 1:3, display = "species", scaling = "none")
  
  # Calcular distancias al centroide (método de Forester et al.)
  dist_centroid <- sqrt(rowSums(loadings^2))
  threshold <- quantile(dist_centroid, probs = 1-alpha)
  
  # Devolver SNPs significativos
  return(names(which(dist_centroid > threshold)))
}

rda_candidates <- signif_snps(rda_result)
rda_candidates

datos_gen_pop <- genind2genpop(datos_gen)
allele_frequencies <- makefreq(datos_gen_pop, quiet = TRUE)
frecuencies_rda_candidates <- allele_frequencies[, rda_candidates]
write.table(frecuencies_rda_candidates,file = "frecuencias_rda.txt",sep = ";",dec = ",",row.names = T)

# Variables

# Scores de las variables ambientales (biplot scores)
var_scores <- scores(rda_result, display = "bp", scaling = "species")

# Calcular correlaciones (equivalentes a loadings en PCA)
correlations <- cor(datos_amb_std, scores(rda_result, display = "sites", choices = 1:2))  # Primeros 2 ejes

# Prueba por términos (variables individuales)
anova_terms <- anova(rda_result, by = "terms", permutations = 10000)

X<-read.csv("nuevo_ambientales.env",sep = ",",dec = ".",header = F)

# Extraer p-values
p_values <- anova_terms$`Pr(>F)`[1:ncol(X)]  # Ignora el último valor (residual)

results_table <- data.frame(
  RDA = "RDA",  # Etiqueta fija o eje específico (ej. "RDA1")
  Variable = colnames(X),
  Correlación = round(correlations[, 1], 3),  # Correlación con RDA1
  p_value = round(p_values, 4)
)

# Añadir asteriscos para significancia
results_table$Significancia <- ifelse(
  results_table$p_value < 0.01, "**",
  ifelse(results_table$p_value < 0.05, "*", "")
)

# Ordenar por correlación absoluta
results_table <- results_table[order(abs(results_table$Correlación), decreasing = TRUE), ]

# Mostrar tabla
print(results_table)

# Combinar correlación y significancia en una columna
results_table$Correlación_formateada <- paste0(
  results_table$Correlación,
  results_table$Significancia
)

# Seleccionar columnas finales
final_table <- results_table[, c("RDA", "Variable", "Correlación_formateada", "p_value")]
colnames(final_table) <- c("RDA", "Variable", "Correlación", "p-value")

# Guardar
write.table(final_table, "resultados_rda.txt", sep=";",dec=",",row.names = TRUE)

# --- 1. Verificar la significancia de los ejes (Opcional pero Recomendado) ---
# Esto ayuda a confirmar que RDA1 es un eje significativo.
print("Resultados de la prueba de significancia por eje:")
anova(rda_result, by = "axis", perm.max = 999)
cat("\n")

# --- 2. Extraer los Loadings de los SNPs (Species Scores) ---
snp_loadings <- scores(rda_result, display = "species")

# --- 3. Identificar y Ordenar los SNPs Extremos por RDA1 ---
# (RDA1 generalmente explica la mayor parte de la variación ambiental)

# a) Obtener las cargas del primer eje de RDA
rda1_loadings <- snp_loadings[, "RDA1"]

# b) Convertir a un dataframe para manipulación
loadings_df <- as.data.frame(rda1_loadings)
loadings_df$SNP <- rownames(loadings_df)
colnames(loadings_df) <- c("RDA1_Loadings", "SNP")

# c) Ordenar por el valor absoluto de la carga (los más extremos primero)
loadings_df_sorted <- loadings_df[order(-abs(loadings_df$RDA1_Loadings)), ]

# d) Imprimir los resultados de los 10 SNPs con mayor loading en RDA1
print("Los 10 SNPs con mayor carga (loading) en el eje RDA1:")
head(loadings_df_sorted, 10)

##

anova_axes <- anova(rda_result, by = "axis", perm.max = 999)
print(anova_axes)

# Determina el número de ejes significativos para usar en la correlación.
# Por simplicidad, usaremos los dos primeros ejes (RDA1 y RDA2),
# que generalmente explican la mayor varianza, pero idealmente se usa
# el número de ejes que resultaron significativos en la prueba ANOVA.
num_ejes_significativos <- 2 # Asume RDA1 y RDA2. Ajusta si es necesario.

# --- Paso 3: Extraer Scores y Calcular la Matriz de Correlación (Proyección) ---

# a) Coordenadas de los SNPs (Species Scores / Loadings)
snp_loadings <- scores(rda_result, display = "species")
# b) Coordenadas de las Variables Ambientales (Environmental Scores / Biplot)
env_scores <- scores(rda_result, display = "bp")

# Seleccionar solo los ejes significativos para el cálculo
snp_coords <- snp_loadings[, 1:num_ejes_significativos, drop = FALSE]
env_coords <- env_scores[, 1:num_ejes_significativos, drop = FALSE]

# Calcular la matriz de correlación/proyección (Producto Escalar)
# Esta matriz mide la fuerza de la asociación de cada SNP con cada variable ambiental
# en el espacio definido por los ejes RDA seleccionados.
correlation_matrix <- snp_coords %*% t(env_coords)

print("\n--- Matriz de Correlación SNP vs. Variables Ambientales ---")
# Muestra las primeras filas y columnas de la matriz resultante
print(head(correlation_matrix))
# Muestra las dimensiones: Filas = SNPs, Columnas = Variables Ambientales
print(paste("Dimensiones de la matriz:", nrow(correlation_matrix), "filas (SNPs) x", ncol(correlation_matrix), "columnas (Variables)"))

# --- Paso 4: Identificar el SNP Más Asociado a Cada Variable ---
print("\n--- SNP Más Fuertemente Asociado a Cada Variable Ambiental (Máxima |Correlación|) ---")
results_list <- list()

for (env_var in colnames(correlation_matrix)) {
  # Obtener los valores de correlación para esta variable
  correlations <- correlation_matrix[, env_var]
  
  # Encontrar el SNP con el valor absoluto más alto (la asociación más fuerte)
  max_cor_index <- which.max(abs(correlations))
  snp_nombre <- rownames(correlation_matrix)[max_cor_index]
  cor_valor <- correlations[max_cor_index]
  
  results_list[[env_var]] <- data.frame(
    Variable = env_var,
    SNP_Asociado = snp_nombre,
    Correlacion_RDA = cor_valor
  )
}

# Combinar los resultados en un solo dataframe
most_associated_snps <- do.call(rbind, results_list)

# Ordenar los resultados por el valor absoluto de la correlación de forma descendente
most_associated_snps_sorted <- most_associated_snps[order(-abs(most_associated_snps$Correlacion_RDA)), ]

# Imprimir el listado final
print(most_associated_snps_sorted)
View(most_associated_snps_sorted)

#### PCA 

amb_data2<-read.csv("nuevo_ambientales.env",sep = ",",dec = ".",header = F)
amb_data2<-amb_data2[-c(1:4),]

pop_label<-stringr::str_sub(rownames(datos_gen@tab),start = 1, end = 4)
pop_label[pop_label=="ICCL"]<-"PPCL"
pop_label[pop_label=="BHAR"]<-"BBAR"
pop_label[pop_label=="BSAR"]<-"BAAR"

ambpop_data <- cbind(amb_data2, pop_label)
colnames(ambpop_data)<-c(names(datos_amb),"pop")
ambpop_data$pop<-as.factor(ambpop_data$pop)

ambpop_data$pop <- factor(ambpop_data$pop,
                         levels=c("LPUY","PEUY","PRUY","BAAR","BBAR","MPAR","FKUK","BACL","PMAR",
                                  "CTCL","RMCL","PCCL","IPCL",
                                  "SGCL","BLCL","PACL","PPCL","BICL","CMCL",
                                  "CAAR","ALAR","PNCL","PWCL"))

#pairs(ambpop_data,lower.panel = NULL,col = as.numeric(ambpop_data$pop))

#Grafico PCA

uru_pca<-prcomp(ambpop_data[,-18],scale. = T)
uru_pca
biplot(uru_pca)

# Coordenadas de sitios (muestras)
sites <- as.data.frame(scores(uru_pca, display = "sites", choices = c(1, 2)))
sites$Species <- ambpop_data$pop  # Añadir columna de grupos

# Coordenadas de especies (variables)
species <- as.data.frame(scores(uru_pca, display = "species", choices = c(1, 2))*10)
species$Variable <- rownames(species)  # Nombres de las variables

# Porcentaje de varianza explicada
var_exp <- round(100 * uru_rda$CA$eig / sum(uru_rda$CA$eig), 1)

#plot
library(ggrepel)

cols = c(colorRampPalette(c("#38F4E5","#64D5F7"))(3),
         colorRampPalette(c("#208A8E","black"))(6),
         colorRampPalette(c("#EC7B74","#BF635A"))(4),
         colorRampPalette(c("#5E8B48","#DAC753"))(10))

ggplot() +
  # Sitios (muestras) coloreados por pop
  geom_point(
    data = sites,
    aes(x = PC1, y = PC2, color = Species),
    size = 3
  ) +
  # Flechas de variables (pop)
  geom_segment(
    data = species,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "black"
  ) +
  # Etiquetas de variables
  geom_text_repel(
    data = species,
    aes(x = PC1, y = PC2, label = Variable),
    color = "black",
    box.padding = 0.8
  ) +
  # Ajustes estéticos
  labs(
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)"),
    #title = ""
  ) +
  scale_color_manual(values = cols) +  # Colores para pop
  theme_minimal() +
  coord_fixed()  # Relación 1:1 entre ejes

ggsave("pca_ambiental.tiff",width = 18,height = 8,dpi = "print")

# Grafico RDA

uru_rda<-rda(ambpop_data[,-18],scale=T)
biplot(uru_rda,display = c("sites","species"),type = c("text","points"))

# Coordenadas de sitios (muestras)
sites <- as.data.frame(scores(uru_rda, display = "sites", choices = c(1, 2)))
sites$Species <- ambpop_data$pop  # Añadir columna de grupos

# Coordenadas de especies (variables)
species <- as.data.frame(scores(uru_rda, display = "species", choices = c(1, 2)))
species$Variable <- rownames(species)  # Nombres de las variables

# Porcentaje de varianza explicada
var_exp <- round(100 * uru_rda$CA$eig / sum(uru_rda$CA$eig), 1)

#plot
library(ggrepel)

cols = c(colorRampPalette(c("#38F4E5","#64D5F7"))(3),
         colorRampPalette(c("#208A8E","black"))(6),
         colorRampPalette(c("#EC7B74","#BF635A"))(4),
         colorRampPalette(c("#5E8B48","#DAC753"))(10))

ggplot() +
  # Sitios (muestras) coloreados por pop
  geom_point(
    data = sites,
    aes(x = PC1, y = PC2, color = Species),
    size = 3
  ) +
  # Flechas de variables (pop)
  geom_segment(
    data = species,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "black"
  ) +
  # Etiquetas de variables
  geom_text_repel(
    data = species,
    aes(x = PC1, y = PC2, label = Variable),
    color = "black",
    box.padding = 0.8
  ) +
  # Ajustes estéticos
  labs(
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)"),
    #title = ""
  ) +
  scale_color_manual(values = cols) +  # Colores para pop
  theme_minimal() +
  coord_fixed()  # Relación 1:1 entre ejes

ggsave("rda_ambiental.tiff",width = 18,height = 8,dpi = "print")

#### Hierfstat

# 2. Análisis FST (requiere paquete hierfstat o similar)
library(hierfstat)
fst <- varcomp.glob(levels = datos_gen$pop, loci = datos_gen$tab)$FST
fst_outliers <- which(fst > quantile(fst, 0.90)) # Percentil 95

# 3. Loci candidatos finales (comunes a los tres métodos)
adaptive_loci <- intersect(rda_candidates, fst_outliers)

#### LFMM

Y<-read.csv("uruprocesado80v2_allref_id_pop.str.lfmm_imputed.lfmm",sep = " ",header=F)
Y<-read.csv("uruprocesado80v2_allref_id_pop.str.lfmm",sep = " ",header=F)
which(colSums(Y == 9) > 0)

X<-read.csv("nuevo_ambientales.env",sep = ",",dec = ".",header=F)
# diferentes

n_filas <- nrow(X)
secuencia_a_sumar <- seq(0.00011, by = 0.00011, length.out = n_filas)

# Sumar la secuencia a todas las columnas numéricas
X <- X %>%
  mutate(across(where(is.numeric), function(x) x + secuencia_a_sumar))

X <- X[,-12]
colnames(X)<-c("Average.OceanTemperature","Average.Salinity","Average.Nitrate","Average.Phosphate","Average.Silicate","Average.SeaIceThickness","Average.SeaIceCover","Average.DissolvedMolecularOxygen","Average.pH","Average.DissolvedIron","Average.Chlorophyll","Average.TotalCloudFraction","Average.AirTemperature")

mod.lfmm2 <- lfmm2(input = Y[-c(1:4),], env = X[-c(1:4),], K = 4)

z_scores <- lfmm2.test(object = mod.lfmm2, input = Y[-c(1:4),], env = X[-c(1:4),])

# Convertir la matriz de p-values a un data.frame largo
results_all <- as.data.frame(z_scores$pvalues) %>%
  rownames_to_column("Variable_ambiental") %>%
  pivot_longer(
    cols = -Variable_ambiental,
    names_to = "SNP",
    values_to = "pvalue"
  ) %>%
  group_by(Variable_ambiental) %>%  # Agrupar por variable ambiental
  mutate(
    # Ajustar p-values por FDR dentro de cada variable
    pvalue_ajustado = p.adjust(pvalue, method = "fdr"),
    # Etiquetar significancia
    Significativo = ifelse(pvalue_ajustado < 0.05, "Sí", "No")
  ) %>%
  ungroup()

# Ejemplo: Vector de nombres reales (ajusta según tu caso)
nombres_reales <- locNames(datos_gen)  # Hasta completar todos tus SNPs

# Verificar que la longitud coincida
stopifnot(length(nombres_reales) == length(unique(results_all$SNP)))

# Crear diccionario (named vector)
diccionario_snps <- setNames(nombres_reales, paste0("V", 1:length(nombres_reales)))


# SNPs significativos para CADA variable ambiental
candidates_by_var <- results_all %>%
  filter(Significativo == "Sí") %>%
  arrange(Variable_ambiental, pvalue_ajustado)

# Ver el resumen
candidates_by_var %>%
  count(Variable_ambiental, name = "N_SNPs_significativos")

print(candidates_by_var, n=Inf)

# Guardar resultados
write_csv(candidates_by_var, "snps_significativos_por_variable.csv")

unique_candidates <- candidates_by_var %>%
  distinct(SNP, .keep_all = TRUE) %>%
  arrange(pvalue_ajustado)