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

#boxplot
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


#Procrustes Alignment
proc_data <- 
  gpagen(epi_data) 


#plot tangent space
pca_data <- 
  plotTangentSpace(proc_data$coords, groups=colour_ids, label=TRUE)


#build geomorph dataframe
gdf <- 
  geomorph.data.frame(
    coords = two.d.array(proc_data$coords), 
    groups = group_ids,
    colors = colour_ids,
    new_stage = factor(curv_data$new_stage),
    Csize = proc_data$Csize)


#is shape determined by developmental stage and taxon?
fit1 <- lm.rrpp(coords ~ new_stage*groups, data = gdf, iter = 199)


#trajectory analysis
TA1 <- 
  trajectory.analysis(
    fit1, 
    groups = gdf$groups, 
    traj.pts = gdf$new_stage, 
    print.progress = FALSE)


#plot groups*stage morphospace
TP1 <- 
  plot(TA1, pch = as.numeric(gdf$groups) + 20, 
       bg = as.numeric(gdf$new_stage),
       cex = 1.5, col = "gray")

#add trajectory lines per species
add.trajectories(TP1, traj.pch = c(21,22), start.bg = 1, end.bg = 2)

#add legend
legend("topright", levels(gdf$groups), pch =  c(21, 22), pt.bg = 1)


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

legend("topright", levels(gdf$groups), col=c(1,2), pch=16, pt.bg = 1)
