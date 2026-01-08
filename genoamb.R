# 
# Description:
#
# Script to perform environmental analysis and test the relationship between ambient and SNP genotypes
#

# Clean Rstudio environment
rm(list=ls()) 
graphics.off() 

#Import libraries
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

#### clone coordinates and environmental data

# coordinates

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

# Create a list to store duplicate rows
filas_duplicadas <- list()

# Iterate through each row of the dataframe and the corresponding number in the vector
for (i in 1:nrow(coordenadas)) {
  fila_a_duplicar <- coordenadas[i, ]  # Select the current row
  num_repeticiones <- conteo$ind[i] # Get the number of times to copy
  
  # Duplicate the row and add it to the list
  for (j in 1:num_repeticiones) {
    filas_duplicadas[[length(filas_duplicadas) + 1]] <- fila_a_duplicar
  }
}

# Convert the list of duplicate rows to a new dataframe
nuevo_coordenadas <- do.call(rbind, filas_duplicadas)
write.table(nuevo_coordenadas,file = "nuevo_coordenadas.txt",sep = ",",dec = ".",
            row.names = F,col.names = F,quote = F)

nuevo_coordenadas2<-as.matrix(nuevo_coordenadas)
nuevo_coordenadas2[, c(1, 2)] <- nuevo_coordenadas2[, c(2, 1)]

plot(nuevo_coordenadas2, pch = 19, cex = .5, 
     xlab = "Longitude (°E)", ylab = "Latitude (°N)")
maps::map(add = T, interior = T)

## Environmental data

ambientales<-read.csv("data_ambiental.env",sep = ",",dec = ".",header = T)

# Verify that the number of rows in the dataframe matches the length of the vector
if (nrow(ambientales) != length(conteo$ind)) {
  stop("The number of rows in the dataframe must be equal to the length of the copy vector.")
}

# Create a list to store duplicate rows
filas_duplicadas <- list()

# Iterate through each row of the dataframe and the corresponding number in the vector
for (i in 1:nrow(ambientales)) {
  fila_a_duplicar <- ambientales[i, ] # Select the current row
  num_repeticiones <- conteo$ind[i] # Get the number of times to copy
  
  # Duplicate the row and add it to the list
  for (j in 1:num_repeticiones) {
    filas_duplicadas[[length(filas_duplicadas) + 1]] <- fila_a_duplicar
  }
}

# Convert the list of duplicate rows to a new dataframe
nuevo_ambientales <- do.call(rbind, filas_duplicadas)
write.table(nuevo_ambientales,file = "nuevo_ambientales.env",sep = ",",dec = ".",
            row.names = F,col.names = F,quote = F)

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

#---- Vegan RDA ----

# Import genotypes
datos_gen <- read.genepop("uruprocesado80v2_onlystudy_nopbuy.gen", ncode = 3)

##### Identify monomorphic SNPs
isPoly(datos_gen) %>% summary

# Exclude monomorphic SNPs
poly_loci = names(which(isPoly(datos_gen) == TRUE))
datos_gen = datos_gen[loc = poly_loci]
isPoly(datos_gen) %>% summary

# Convert to allelic frequencies
frec_alelicas <- genind2genpop(datos_gen) %>% tab()

# Upload environmental data (must have the same samples as the genetic data)
datos_amb <- read.csv("data_ambiental.txt", row.names = 1)
nombres_nuevos <- sub("^[^.]+\\.", "", names(datos_amb))
names(datos_amb) <- nombres_nuevos
datos_amb <- datos_amb[-c(3,19),]

rownames(frec_alelicas)<-rownames(datos_amb)

# Make sure the rows match
frec_alelicas <- frec_alelicas[rownames(datos_amb), ]

# Standardize environmental data (important for RDA)
datos_amb_std <- scale(datos_amb)
datos_amb_std<-as.data.frame(datos_amb_std)

# Delete columns with NA

datos_amb_std2 <- datos_amb_std %>%
  select(where(~ !anyNA(.)))

# Correlation between environmental variables

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

# Calculate full correlation matrix (SNPs × variables)
cor_matrix <- cor(frec_alelicas, datos_amb_std2,method = "spearman")

library(Hmisc)

datos_combinados <- cbind(frec_alelicas, datos_amb_std2)
resultados_cor <- rcorr(as.matrix(datos_combinados), type = "spearman")
p_valores <- resultados_cor$P
nombres_marcadores <- colnames(frec_alelicas)
nombres_ambientales <- colnames(datos_amb_std2)

p_valores_filtrados <- p_valores[nombres_marcadores, nombres_ambientales]
print(p_valores_filtrados)

# Obtain the filtered p-value matrix (as we did before)
library(Hmisc)
library(dplyr)
library(tidyr)
library(tibble)

datos_combinados <- cbind(frec_alelicas, datos_amb_std2)
resultados_cor <- rcorr(as.matrix(datos_combinados), type = "spearman")
p_valores_filtrados <- resultados_cor$P[colnames(frec_alelicas), colnames(datos_amb_std2)]

# Transform the p-value matrix to long format
p_valores_largos <- p_valores_filtrados %>%
  as.data.frame() %>%
  rownames_to_column("SNP") %>%
  pivot_longer(
    cols = -SNP,
    names_to = "Predictor",
    values_to = "p_value"
  )

# Reorganize the correlation matrix into long format and combine it with the p-values
resultados_finales <- resultados_cor$r[colnames(frec_alelicas), colnames(datos_amb_std2)] %>%
  as.data.frame() %>%
  rownames_to_column("SNP") %>%
  pivot_longer(
    cols = -SNP,
    names_to = "Predictor",
    values_to = "Correlation"
  ) %>%
  # Join with the p-value table using SNP and Predictor as keys
  left_join(p_valores_largos, by = c("SNP", "Predictor")) %>%
  mutate(
    N = 1:n()
  ) %>%
  # Delete the Eje_RDA column
  select(N, SNP, Predictor, Correlation, p_value) %>%
  arrange(Predictor, desc(abs(Correlation))) %>%
  mutate(
    Correlation = round(Correlation, 2),
    p_value = round(p_value, 4) # Round p-values
  )

print(resultados_finales)

# Filter correlations with |r| > 0.8 (adjust this threshold)
cor_extremas <- resultados_finales %>%
  filter(abs(p_value) < 0.05 & abs(Correlation) >= 0.8) %>%
  arrange(desc(abs(p_value)))

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

# Identify significant SNPs (p < 0.05)
signif_snps <- function(rda_obj, alpha = 0.05){
  loadings <- scores(rda_obj, choices = 1:3, display = "species", scaling = "none")
  dist_centroid <- sqrt(rowSums(loadings^2))
  threshold <- quantile(dist_centroid, probs = 1-alpha)
  return(names(which(dist_centroid > threshold)))
}

rda_candidates <- signif_snps(rda_result)
rda_candidates

datos_gen_pop <- genind2genpop(datos_gen)
allele_frequencies <- makefreq(datos_gen_pop, quiet = TRUE)
frecuencies_rda_candidates <- allele_frequencies[, rda_candidates]
write.table(frecuencies_rda_candidates,file = "frecuencias_rda.txt",sep = ";",dec = ",",row.names = T)

# Variables

# Environmental variable scores (biplot scores)
var_scores <- scores(rda_result, display = "bp", scaling = "species")

# Calculate correlations (equivalent to loadings in PCA)
correlations <- cor(datos_amb_std, scores(rda_result, display = "sites", choices = 1:2))  # Primeros 2 ejes

# Test by terms (individual variables)
anova_terms <- anova(rda_result, by = "terms", permutations = 10000)

X<-read.csv("nuevo_ambientales.env",sep = ",",dec = ".",header = F)

# Extract p-values
p_values <- anova_terms$`Pr(>F)`[1:ncol(X)]

results_table <- data.frame(
  RDA = "RDA",  
  Variable = colnames(X),
  Correlación = round(correlations[, 1], 3), 
  p_value = round(p_values, 4)
)

# Add asterisks for significance
results_table$Significancia <- ifelse(
  results_table$p_value < 0.01, "**",
  ifelse(results_table$p_value < 0.05, "*", "")
)

# Order by absolute correlation
results_table <- results_table[order(abs(results_table$Correlación), decreasing = TRUE), ]
print(results_table)

# Combine correlation and significance in one column
results_table$Correlación_formateada <- paste0(
  results_table$Correlación,
  results_table$Significancia
)

# Select final columns
final_table <- results_table[, c("RDA", "Variable", "Correlación_formateada", "p_value")]
colnames(final_table) <- c("RDA", "Variable", "Correlación", "p-value")

# Save table
write.table(final_table, "resultados_rda.txt", sep=";",dec=",",row.names = TRUE)

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

# PCA plot
uru_pca<-prcomp(ambpop_data[,-18],scale. = T)
uru_pca
biplot(uru_pca)

# Site coordinates (samples)
sites <- as.data.frame(scores(uru_pca, display = "sites", choices = c(1, 2)))
sites$Species <- ambpop_data$pop 

# Species coordinates (variables)
species <- as.data.frame(scores(uru_pca, display = "species", choices = c(1, 2))*10)
species$Variable <- rownames(species)  
                
# Percentage of variance explained
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

# RDA plot
uru_rda<-rda(ambpop_data[,-18],scale=T)
biplot(uru_rda,display = c("sites","species"),type = c("text","points"))

# Site coordinates (samples)
sites <- as.data.frame(scores(uru_rda, display = "sites", choices = c(1, 2)))
sites$Species <- ambpop_data$pop  # Añadir columna de grupos

# Species coordinates (variables)
species <- as.data.frame(scores(uru_rda, display = "species", choices = c(1, 2)))
species$Variable <- rownames(species)  # Nombres de las variables

# Percentage of variance explained
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
