# 
# Description:
#
# Script to clean genotypes by SNPs and individuals, obtein different genetic stadistics (FIS, FST, Ho and He values) and to perform the DAPC analysis using the library Adegenet 
#

# Clean Rstudio environment
rm(list=ls()) 
graphics.off() 

# Import libraries
library(tidyverse)
library(adegenet)
library(poppr)
library(ggrepel)
library(hierfstat)

setwd("D:/Trabajo/Paper/uruguay/adegenet/")
setwd("/media/oxus/LS/Trabajo/Paper/uruguay/adegenet")

datos_ori<-read.genepop("uruguay.gen",ncode = 3,quiet=FALSE)
datos_ori<-read.genepop("uruprocesado80v2.gen",ncode = 3,quiet = FALSE)

# Identify SNPs removed during import

nLoc(datos_ori)
snp<-read.table(file="clipboard",sep = "\t",header = F) #Original SNPs list from genind file
setdiff(snp$V1,locNames(datos_ori)) #SNP not imported
setdiff(locNames(dato_ori),snp$V1)

# Edit population names
unique(datos_ori$pop)
popNames(datos_ori)<-c("PBUY","LPUY","PEUY","PRUY","BSAR","BBAR","MPAR","FKUK","BACL","PMAR",
                       "CTCL","RMCL","PCCL","IPCL",
                       "SGCL","BLCL","PACL","PPCL","BICL","CMCL",
                       "CAAR","ALAR","PNCL","PWCL",
                       "Mch","Mp","Mg_Med","Mg_Atl","Me_AM","Me_EU","Mpl","Mt")
unique(datos_ori$pop)

#### Save a copy of the original dataset without modification

datos = datos_ori

#### Basic Stadistics (N° per population, N° alelles, etc etc)

summary(datos_ori$pop)
summary(datos_ori)

#### Basic Stadistics with Hierfstat

hierfstat_data <- genind2hierfstat(datos)

# Using adegenet::summary() for a quick overview
adegenet_summary <- summary(datos)

n_individuos_por_poblacion <- as.data.frame(table(pop(datos_study))) %>%
  rename(Poblacion = Var1, n = Freq)

stats_result <- basic.stats(hierfstat_data_study)

#Obtein Ho, He and FIS values per population
ho_promedio_poblacion <- colMeans(stats_result$Ho, na.rm = TRUE)
hs_promedio_poblacion <- colMeans(stats_result$Hs, na.rm = TRUE)
fis_promedio_poblacion <- colMeans(stats_result$Fis, na.rm = TRUE)

# Generate table with the previously obtained values
tabla_resumen_poblacion <- tibble(
  Poblacion = names(ho_promedio_poblacion),
  Ho = ho_promedio_poblacion,
  Hs = hs_promedio_poblacion,
  Fis = fis_promedio_poblacion,
)

# Add the n per population 
tabla_resumen_poblacion_final <- tabla_resumen_poblacion %>%
  left_join(n_individuos_por_poblacion, by = "Poblacion") %>%
  select(Poblacion, n, everything())

# Made a results table
write.table(tabla_resumen_poblacion_final,file = "tabla_resumen_poblacion_final_study.txt", sep = ",", dec = ".",row.names = F)
write.table(stats_result$overall,file = "stats_result_study.txt", sep = ",", dec = ".",row.names = F)

# Exclude SNPs and samples with a genotyping rate below 80%

## by loci

locmiss_datos = propTyped(datos, by = "loc")
mean(locmiss_datos)

locmiss_datos[which(locmiss_datos < 0.80)] # print loci with < XX% complete genotypes ## named numeric(0)
mean(locmiss_datos[which(locmiss_datos < 0.80)])
missingno(datos, type = "loci", cutoff = 0.20)
#barplot(locmiss_datos, ylim = c(0,1), ylab = "Complete genotypes (proportion)", xlab = "Locus", las = 2, cex.names = 0.5)

datos = missingno(datos, type = "loci", cutoff = 0.20)

### by individuals

indmiss_datos = propTyped(datos, by = "ind")
mean(indmiss_datos)
indmiss_datos[ which(indmiss_datos < 0.80) ] # print individuals with < 70% complete genotypes 
mean(indmiss_datos[ which(indmiss_datos < 0.80) ])
datos80 = missingno(datos, type = "geno", cutoff = 0.20)

### Identify monomorphic SNPs
isPoly(datos) %>% summary

datos_study = popsub(datos, sublist = c("PBUY","LPUY","PEUY","PRUY","BSAR","BBAR","MPAR","FKUK","BACL","PMAR","CTCL",
                                        "RMCL","PCCL","IPCL","SGCL","BLCL","PACL","PPCL","BICL","CMCL","CAAR","ALAR",
                                        "PNCL","PWCL")) 
isPoly(datos_study) %>% summary

### Exclude monomorphic SNPs

poly_loci = names(which(isPoly(datos_study) == TRUE))
datos_study = datos_study[loc = poly_loci]
isPoly(datos_study) %>% summary

hierfstat_data_study <- genind2hierfstat(datos_study)
basic.stats(hierfstat_data_study)

# Create a genepop file with the processed genetic data

devtools::source_url("https://raw.githubusercontent.com/Tom-Jenkins/utility_scripts/master/TJ_genind2genepop_function.R")
genind2genepop(datos80,file = "uruprocesado80v21.gen")

######
summary(datos80$pop)
as.data.frame(table(datos80$pop))

##### Obtein FST matrix and create a graph of this
datos80=datos
datos_fst = genet.dist(datos80, method = "WC84") %>% round(digits = 4)
write.table(as.matrix(datos_fst),file = "fst.txt", sep = ",", dec = ".",row.names = F)

# plot FST values
lab_order = c(unique(datos80$pop)) # Ordenar poblaciones

# Change order of rows and cols
fst.mat = as.matrix(datos_fst)
fst.mat1 = fst.mat[lab_order, ]
fst.mat2 = fst.mat1[, lab_order]

# Create a data.frame
ind = which(upper.tri(fst.mat2), arr.ind = TRUE)
fst.df = data.frame(Site1 = dimnames(fst.mat2)[[2]][ind[,2]],
                    Site2 = dimnames(fst.mat2)[[1]][ind[,1]],
                    Fst = fst.mat2[ ind ])

# Keep the order of the levels in the data.frame for plotting 
fst.df$Site1 = factor(fst.df$Site1, levels = unique(fst.df$Site1))
fst.df$Site2 = factor(fst.df$Site2, levels = unique(fst.df$Site2))

fst.df$Fst[fst.df$Fst < 0] = 0 # Convert minus values to zero
fst.df %>% str # Print data.frame summary
fst.label = expression(italic("F")[ST]) # Fst italic label
mid = max(fst.df$Fst) / 2 # Extract middle Fst value for gradient argument

# Plot FST
ggplot(data = fst.df, aes(x = Site1, y = Site2, fill = Fst)) +
  geom_tile(colour = "black") +
  geom_text(aes(label = Fst), color="black", size = 3) +
  scale_fill_gradient2(low = "#74ADD1", mid = "grey90", high = "#F46D43", midpoint = mid, name = fst.label, limits = c(0, max(fst.df$Fst)), breaks = c(0, 0.5, 1)) +
  scale_x_discrete(expand = c(0,0),limits=rev) +
  scale_y_discrete(expand = c(0,0), position = "left") +
  theme(axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5, face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        #legend.title = element_text(size = 14, face = "bold"),
        #legend.text = element_text(size = 10)
  )

ggsave("fst.tiff",width = 18,height = 8,dpi = "print")

####### Prepare different groups to analyze

datos80 = datos_ori

#### Only references

datos = popsub(datos80, sublist = c("Mch","Mp","Mg_Med","Mg_Atl","Me_AM","Me_EU","Mpl","Mt"))

#### Only studied locations

datos = popsub(datos80, exclude = c("Mch","Mp","Mg_Med","Mg_Atl","Me_AM","Me_EU","Mpl","Mt"))

#### mch and mp references + studied locations

datos = popsub(datos80, exclude = c("Mg_Med","Mg_Atl","Me_AM","Me_EU","Mpl","Mt"))

#### mch, mg y mp + studied locations

datos = popsub(datos80, exclude = c("Mg_Atl","Me_AM","Me_EU","Mpl","Mt"))

###### DAPC analysis

set.seed(123)
x = tab(datos, NA.method = "mean")
crossval = xvalDapc(x, datos$pop, result = "groupMean", xval.plot = TRUE)

crossval$`Root Mean Squared Error by Number of PCs of PCA` # Number of PCs with best stats (lower score = better)
crossval$`Number of PCs Achieving Highest Mean Success`
crossval$`Number of PCs Achieving Lowest MSE`
numPCs = as.numeric(crossval$`Number of PCs Achieving Lowest MSE`)

grup_referenciayambos <- find.clusters(datos, n.pca = numPCs, max.n.clust=30) # Run a DAPC using a specific k number
dapc1 = dapc(datos, grup_referenciayambos$grp, n.pca = nPop(datos)-1, n.da = length(levels(grup_referenciayambos$grp))-1)

scatter(dapc1)
table(pop(datos),dapc1$assign)

#### Plot by population

# Analyse how much percent of genetic variance is explained by each axis
percent = dapc1$eig/sum(dapc1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,80),
        names.arg = round(percent, 1),las=1)

ind_coords = as.data.frame(dapc1$ind.coord) # Create a data.frame containing individual coordinates

colnames(ind_coords) = c("Axis1","Axis2") # Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2","Axis3") # Rename columns of dataframe

ind_coords$Ind = indNames(datos) # Add a column containing individuals
ind_coords$Site = datos$pop # Add a column with the site IDs

centroid = aggregate(cbind(Axis1, Axis2) ~ Site, data = ind_coords, FUN = mean) # Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords, FUN = mean) # Calculate centroid (average) position for each population

ind_coords = left_join(ind_coords[,-c(4,5)], centroid, by = "Site", suffix = c("",".cen")) # Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen")) # Add centroid coordinates to ind_coords dataframe

# Colors by different groups

# Only references

cols = c("indianred1","plum","goldenrod1","goldenrod3","dodgerblue","dodgerblue4",
         "darkolivegreen4","grey60")

cols = c("indianred1","plum","goldenrod1","dodgerblue","darkolivegreen4","grey60")

# Only studied locations

cols = c(colorRampPalette(c("#38F4E5","#64D5F7"))(4),
         colorRampPalette(c("#208A8E","#306478"))(6),
         colorRampPalette(c("#EC7B74","#BF635A"))(4),
         colorRampPalette(c("#5E8B48","#DAC753"))(10))

# 3 references + studied locations

cols = c(colorRampPalette(c("#38F4E5","#64D5F7"))(4),
         colorRampPalette(c("#208A8E","#306478"))(6),
         colorRampPalette(c("#EC7B74","#BF635A"))(4),
         colorRampPalette(c("#5E8B48","#DAC753"))(10),
         "indianred1","plum","goldenrod1")

# 5 references + studied locations

cols = c(colorRampPalette(c("#38F4E5","#64D5F7"))(4),
         colorRampPalette(c("#208A8E","#306478"))(6),
         colorRampPalette(c("#EC7B74","#BF635A"))(4),
         colorRampPalette(c("#5E8B48","#DAC753"))(10),
         "indianred1","plum","goldenrod1","dodgerblue","darkolivegreen4")

# 6 references + studied locations

cols = c(colorRampPalette(c("#38F4E5","#64D5F7"))(4),
         colorRampPalette(c("#208A8E","#306478"))(6),
         colorRampPalette(c("#EC7B74","#BF635A"))(4),
         colorRampPalette(c("#5E8B48","#DAC753"))(10),
         "indianred1","plum","goldenrod1","dodgerblue","darkolivegreen4","grey60")

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Scatter plot axis 1 vs. 2 solo lineas
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  #geom_point(aes(color = Site), shape = c(rep(1,4),rep(2,29),rep(0,24),rep(5,20),
   #                                       rep(16,20),rep(17,20),rep(15,20),rep(18,13),
    #                                      rep(21,19),rep(22,17),rep(23,27),rep(24,19),rep(25,20),rep(8,28),
     #                                     rep(6,27),rep(3,23),rep(4,20),rep(1,29),
      #                                    rep(16,30),rep(17,30),rep(15,36),rep(18,20),rep(19,21),rep(7,20)
       #                                   ),
        #     size = 3, show.legend = FALSE)+
  # centroids con etiquetas separadas
  geom_label_repel(data = centroid, aes(label = Site, fill = Site), size = 3, show.legend = FALSE,max.overlaps = 100)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  ggtitle("")+
  #coord_fixed()+
  theme(axis.text.y = element_text(colour="black", size=12),
        axis.text.x = element_text(colour="black", size=12),
        axis.title = element_text(colour="black", size=12),
        panel.border = element_rect(colour="black", fill=NA),
        panel.background = element_blank(),
        plot.title = element_text(hjust=0.5, size=15)
  )

# Scatter plot axis 1 vs. 2 solo puntos
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  #geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  geom_point(aes(color = Site), shape = c(rep(15,4),rep(15,29),rep(15,24),rep(15,20), #uru
                                         rep(16,30),rep(16,30),rep(16,36), #arg
                                         rep(17,20),rep(17,21),rep(17,20), #fkuk bacl pmar
                                         rep(18,27),rep(18,23),rep(18,20),rep(18,29), # ctcl rmcl pccl ipcl
                                         rep(17,20),rep(17,20),rep(17,20),rep(17,13), #magallanes
                                         rep(17,19),rep(17,17),rep(17,27),rep(17,19),rep(17,20),
                                         rep(17,28)
                                     ),
       size = 3, show.legend = FALSE)+
  # centroids con etiquetas separadas
  geom_label_repel(data = centroid, aes(label = Site, fill = Site), size = 3, show.legend = FALSE,max.overlaps = 100)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  ggtitle("")+
  #coord_fixed()+
  theme(axis.text.y = element_text(colour="black", size=12),
        axis.text.x = element_text(colour="black", size=12),
        axis.title = element_text(colour="black", size=12),
        panel.border = element_rect(colour="black", fill=NA),
        panel.background = element_blank(),
        plot.title = element_text(hjust=0.5, size=15)
  )


# Scatter plot axis 1 vs. 2 solo puntos con elipse

k_groups <- dapc1$grp
ind_coords$K_Group <- k_groups
centroid_k <- aggregate(ind_coords[, c("Axis1", "Axis2")], by = list(K_Group = ind_coords$K_Group), FUN = mean)

ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  stat_ellipse(
    aes(color = K_Group, group = K_Group), # <-- CAMBIO CLAVE: Usar K_Group
    type = "t", 
    level = 0.95, 
    linewidth = 1,
    show.legend = TRUE # Mostrar leyenda de grupos K
  ) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  #geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  geom_point(aes(color = Site), shape = c(rep(15,4),rep(15,29),rep(15,24),rep(15,20), #uru
                                          rep(16,30),rep(16,30),rep(16,36), #arg
                                          rep(17,20),rep(17,21),rep(17,20), #fkuk bacl pmar
                                          rep(18,27),rep(18,23),rep(18,20),rep(18,29), # ctcl rmcl pccl ipcl
                                          rep(17,20),rep(17,20),rep(17,20),rep(17,13), #magallanes
                                          rep(17,19),rep(17,17),rep(17,27),rep(17,19),rep(17,20),
                                          rep(17,28)
  ),
  size = 3, show.legend = FALSE)+
  # centroids con etiquetas separadas
  geom_label_repel(
    data = centroid_k, # <-- Usar los centroides de K_Group
    aes(label = K_Group, fill = K_Group), 
    size = 3, 
    show.legend = FALSE,
    max.overlaps = 100
  ) +
  # colouring
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = c("black","black","black",cols)) +
  # custom labels
  labs(x = xlab, y = ylab)+
  ggtitle("")+
  #coord_fixed()+
  theme(axis.text.y = element_text(colour="black", size=12),
        axis.text.x = element_text(colour="black", size=12),
        axis.title = element_text(colour="black", size=12),
        panel.border = element_rect(colour="black", fill=NA),
        panel.background = element_blank(),
        plot.title = element_text(hjust=0.5, size=15)
  )

# Save DAPC plot

ggsave("onlystudyv3_dot.tiff",width = 18,height = 8,dpi = "print")

# Obtener las probabilidades a posteriori


write.table(round(dapc1$posterior,digits = 4), file = "pos_3ref.txt",sep = "\t")
