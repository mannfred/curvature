library(tidyverse)
library(here)
library(geomorph)


#read together
epi_data <- 
  readmulti.tps(
  c(here("data/epimedium_photos/koreanum/koreanum_dorsal_appended_geomorph.TPS"),
    here("data/epimedium_photos/violaceum/violaceum_dorsal_appended_geomorph.TPS")),
  specID="imageID"
)

#import size class data
#rename and bin developmental stages
curv_data <- 
  read.csv(here("data/epimedium_curv_size_data.csv")) %>% 
  rename(sample = species_individual_panicle_flower) %>% 
  slice(20:77) %>%  #keep E. kor and E. vio
  mutate(new_stage = case_when(
    stage == "E" ~ "C",
    stage == "G" ~ "G",
    stage == "ND" | stage == "F" | stage == "N" | stage == "O" ~ "T",
    stage == "P" | stage == "D" | stage == "A" ~ "A"))

#boxplot: stage vs size
ggplot(
  data=curv_data %>% 
    filter(species=="koreanum" | species=="violaceum"), #remove filter to include grandiflorum
  aes(x=new_stage, y=sepal_size_mm)) +
geom_boxplot(
  aes(fill = factor(species), #colour by species
      factor(new_stage, levels=c("C", "G", "T", "A")) #reorder
    )
  ) +
  theme(axis.text.x = element_text(angle=90)) +
  theme_classic() + #removes gray backdrop
  theme(legend.position="bottom")  



#groups
group_ids <-
  c(replicate(31, "E. koreanum"),
    replicate(27, "E. violaceum")) %>% 
  factor()

#colours from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/ 
colour_ids <-
  c(replicate(31, "#E69F00"), #31 E. koreanum samples
    replicate(27, "#56B4E9")) #27 E. violaceum samples

#Procrustes Alignment
proc_data <- 
  gpagen(epi_data) 


#plot tangent space
pca_data <- 
  plotTangentSpace(
    proc_data$coords, 
    groups=group_ids, 
    label=TRUE, 
    legend=TRUE, 
    col=colour_ids)


#build geomorph dataframe

gdf <- 
  geomorph.data.frame(
    coords = two.d.array(proc_data$coords), 
    groups = group_ids,
    colors = colour_ids,
    new_stage = factor(curv_data$new_stage),
    Csize = proc_data$Csize,
    #give stage classes a numerical value 
    #(so that traj analysis doesn't order them alphabetically)
    numerical_stage = case_when(curv_data$new_stage == 'C' ~ 1, 
                                curv_data$new_stage == 'G' ~ 2,
                                curv_data$new_stage == 'T' ~ 3,
                                curv_data$new_stage == 'A' ~ 4))

#to play nice with lm
gdf$numerical_stage <- factor(gdf$numerical_stage) 


saveRDS(gdf, file="geomorph_data_frame.rds")



######################################################
#is shape determined by developmental stage and taxon?
fit1 <- lm.rrpp(coords ~ numerical_stage*groups, data = gdf, iter = 199)


#trajectory analysis
TA1 <- 
  trajectory.analysis(
    fit1, 
    groups = gdf$groups, 
    traj.pts = gdf$numerical_stage, 
    print.progress = FALSE)


#plot groups*stage morphospace

traj_col <- 
  gdf$numerical_stage %>% 
  enframe() %>% 
  mutate(colour = 
           case_when(value == 1 ~ "#009E73", #green
                     value == 2 ~ "#CC79A7", #pink
                     value == 3 ~ "#56B4E9", #blue
                     value == 4 ~ "#E69F00",)) #orange



#plotting trajectory analysis AND warp grids
#divide up the plotting window into a 3×3 grid, 
#the numbers correspond to the order and location of each item being plotted
mat <- matrix(c(4,5,0,1,1,2,1,1,3), 3) 

# set the size of the rows and columns
layout(mat, widths=c(1,1,1), heights=c(1,1,0.6)) 

# sets the margins for traj plot
par(mar=c(4, 4, 1, 1))

TP1 <- 
  plot(TA1, pch = as.numeric(gdf$groups) +20, 
       bg = traj_col$colour,
       cex = 2,
       col='white')

#add trajectory lines per species
add.trajectories(
  TP1, 
  traj.pch = c(1,0), 
  traj.lty=c(1,3)
  )

#add legend
legend("topleft",
       levels(gdf$groups), 
       pch =  c(1,0), 
       pt.bg = 1)

legend("topright",
       levels(gdf$numerical_stage), 
       pch = 18, 
       title = "Stages",
       col = c("#009E73", "#CC79A7", "#56B4E9", "#E69F00"))

#add warp grids
ref <- mshape(proc_data$coords)

par(mar = c(0,0,0,0))
plotRefToTarget(ref, pca_data$pc.shapes$PC1min) #PC1 min
plotRefToTarget(ref, pca_data$pc.shapes$PC1max) #PC1 max
plotRefToTarget(ref, pca_data$pc.shapes$PC2min) #PC2 min
plotRefToTarget(ref, pca_data$pc.shapes$PC1max) #PC2 max




###############################################################
#Does shape differ between taxa at a given developmental stage?
fit3 <- procD.lm(coords ~ new_stage*groups, data = gdf, iter =199, RRPP=TRUE)

epi_groups <- interaction(gdf$groups, gdf$new_stage)

pw1<- pairwise(fit3, groups=epi_groups) 

summary(pw1, confidence=0.95, test.type="dist")




###########
# allometry
fit2 <- procD.lm(coords ~ Csize*groups, data = gdf, iter = 199)

plotAllometry(
  fit2, 
  size=gdf$Csize, 
  logsz = TRUE, 
  method="PredLine",
  col=as.numeric(gdf$groups),
  cex=1.5,
  pch=16)

legend("topright", 
       levels(gdf$groups), 
       col=c(1,2), pch=16, pt.bg = 1)
