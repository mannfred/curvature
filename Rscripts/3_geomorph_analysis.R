library(tidyverse)
library(here)
library(geomorph)

#https://qualityandinnovation.com/2019/09/26/why-anova-and-linear-regression-are-the-same/
#http://www.emmasherratt.com/uploads/2/1/6/0/21606686/quick_guide_to_geomorph-_introduction.html
#https://www.cambridge.org/core/journals/paleobiology/article/many-faces-of-synapsid-cranial-allometry/D4DF84985CA3222FBE327AD993394E3D/core-reader


#read separately 
kor_data <- 
  readland.tps(
    here("data/epimedium_photos/koreanum/koreanum_dorsal_appended_geomorph.TPS"), 
    specID="imageID")

vio_data <- 
  readland.tps(
    here("data/epimedium_photos/violaceum/violaceum_dorsal_appended_geomorph.TPS"), 
    specID="imageID")


#read together
epi_data <- 
  readmulti.tps(
  c(here("data/epimedium_photos/koreanum/koreanum_dorsal_appended_geomorph.TPS"),
    here("data/epimedium_photos/violaceum/violaceum_dorsal_appended_geomorph.TPS")),
  specID="imageID"
)

#colours
colour_ids <-
  c(replicate(31, "#3e8dba"), #31 E. koreanum samples
    replicate(27, "#663593")) #27 E. violaceum samples

#groups
group_ids <-
  c(replicate(31, "E. koreanum"),
    replicate(27, "E. violaceum"))

#Procrustes Alignment
proc_data <- gpagen(epi_data) 

#basic PCA
pca_data <- plotTangentSpace(proc_data$coords, groups=colour_ids, label=TRUE)

#build geomorph dataframe
gdf <- 
  geomorph.data.frame(
    proc_data, 
    PC1 = pca_data$pc.scores[,1], 
    PC2 = pca_data$pc.scores[,2], 
    PC3 = pca_data$pc.scores[,3],
    groups = group_ids, 
    colors=colour_ids)


  

fit1 <- procD.lm(coords ~ log(Csize)*groups, data=gdf, iter=999, RRPP=FALSE)
plotAllometry(
  fit1, 
  size=gdf$Csize, 
  logsz = TRUE, 
  method="PredLine",
  col=colour_ids,
  cex=1.5,
  pch=16,
  alpha = I(1/2))




#PC variances
pca_data$pc.summary


#Test whether species differ in ontogenetic trajectory: All Species
AllTaxa.OntogenyComp <- procD.allometry(coords ~ log(Csize), ~species, logsz = TRUE, data=gdf, iter = 9999)



curv_data <- read.csv(here("data/epimedium_curv_size_data.csv")) %>% 
             rename(sample = species_individual_panicle_flower) %>% 
             filter(grepl('K', sample)) #filter for E. koreanum

proc_data$size <- curv_data$sepal_size_mm

fit <- lm.rrpp(coords ~ Csize, data=proc_data, iter=199)

data(plethodon)
plethodon$site

gpagen(lmk_data) %>%
  procD.lm()








#exponential growth
y <- c(1,2,4,5,7,8,9,10,11,13,14,20,24,27,30,34,47,54,60,67,78,96,117,157,196,250,324,425,569,690,846,1055,1281,1429,2061,2541)
x <- 1:36
data<-data.frame(x,y)
mod <- nls(y ~ exp(a + b * x), data = data, start = list(a = 1, b = 1))
mod

plot(data)
lines(data$x, predict(mod))

f <- function(x) exp(-0.2571 + 0.2246*x)
#f(40) = 6166 cases predicted for March 28
