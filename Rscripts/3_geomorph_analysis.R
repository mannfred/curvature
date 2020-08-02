#if: "Error: maybe not installed for this architecture?" then:
#remotes::install_github('github package directory',
#                        INSTALL_opts = c("--no-multiarch")
#)

library(tidyverse)
library(here)
library(geomorph)
library(lme4)
library(lmerTest)
here=here::here

# --------------------
# data import


# import shape data 
epi_lmk_data <- 
  readmulti.tps(
  c(here("data/epimedium_photos/koreanum/koreanum_appended_geomorph.TPS"),
    here("data/epimedium_photos/violaceum/violaceum_appended_geomorph.TPS")),
  specID="imageID")[,,-31] #remove K_1_2_1 (no curvature data)


# import size class data
# rename and bin developmental stages
sizeclass_data <- 
  read.csv(here("data/epimedium_curv_size_data_nogran.csv")) %>% 
  rename(sample = species_individual_panicle_flower) %>% 
  mutate(new_stage = case_when(
    stage == "E" ~ "C",
    stage == "G" ~ "G",
    stage == "ND" | stage == "F" | stage == "N" | stage == "O" ~ "T",
    stage == "P" | stage == "D" | stage == "A" ~ "A")) 


# check that IDs match
identical(as.character(sizeclass_data$sample), dimnames(epi_lmk_data)[[3]])


# plot stage vs size
ggplot(
  data=sizeclass_data %>% 
    filter(species=="koreanum" | species=="violaceum"), #remove filter to include grandiflorum
  aes(x=new_stage, y=sepal_size_mm)) +
geom_boxplot(
  aes(fill = factor(species), #colour by species
      factor(new_stage, levels=c("C", "G", "T", "A")))) + #reorder
  theme(axis.text.x = element_text(angle=90)) +
  theme_classic() + #removes gray backdrop
  theme(legend.position="bottom")  


# Generalized Procrustes Analysis
# First performs a partial Procrustes Superimposition 
# shapes are scaled to unit-centroid size
proc_data <- 
  gpagen(epi_lmk_data) 

# saveRDS(proc_data, file=here("data/RDS_files/proc_data.rds"))

# illustrates that there is shape variation
plot(proc_data) 

# p.147 of APRIMER
# the 58 sets of landmark coordinates need to be 
# converted to 'shape variables', and these shape variables
# are subjected to PCA
# plotTangentSpace() first converts data using two.d.array(),
# then performs prcomp() and plot()

# groups
group_ids <-
  c(replicate(30, "E. koreanum"),
    replicate(27, "E. violaceum")) %>% 
  factor()

# colours from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/ 
colour_ids <-
  c(replicate(30, "#E69F00"), #31 E. koreanum samples
    replicate(27, "#56B4E9")) #27 E. violaceum samples

# PCA
pca_data <- 
  plotTangentSpace(
    proc_data$coords, 
    groups=group_ids, 
    label=TRUE, 
    legend=TRUE, 
    col=colour_ids)

# saveRDS(pca_data, file=here("data/RDS_files/pca_data.rds"))



# ----------------------------
# build geomorph dataframe for
# trajectory analysis


# build geomorph dataframe
# give stage classes a numerical value 
# (so that traj analysis doesn't order them alphabetically)
gdf <- 
  geomorph.data.frame(
    coords = two.d.array(proc_data$coords),
    indiv  = factor(sizeclass_data$indiv),
    taxon = group_ids,
    colors = colour_ids,
    new_stage = factor(sizeclass_data$new_stage),
    Csize = proc_data$Csize,
    numerical_stage = 
      case_when(
        sizeclass_data$new_stage == 'C' ~ 1, 
        sizeclass_data$new_stage == 'G' ~ 2,
        sizeclass_data$new_stage == 'T' ~ 3,
        sizeclass_data$new_stage == 'A' ~ 4))


#to play nice with lm()
gdf$numerical_stage <- factor(gdf$numerical_stage) 

# save for later use
# saveRDS(gdf, file="geomorph_data_frame.rds")



# -------------------------------------
# modelling shape ~ developmental stage


# is shape determined by developmental stage and taxon?
fit1 <- 
  lm.rrpp(coords ~ numerical_stage*taxon + indiv,
          data = gdf,
          SS.type = "I",
          iter = 999)

# developmental stage and taxon*indiv are
# signficant predictors of shape 
# indiv is a random effect
anova.lm.rrpp(
  fit1, 
  effect.type="F", 
  error = c("Residuals", "Residuals", "indiv", "Residuals"))


# make pairwise comparisons of shape for every developmental stage
epi_groups <- interaction(gdf$numerical_stage, gdf$taxon)

pw1 <- RRPP::pairwise(fit1, groups = epi_groups) 

sum <- RRPP::summary.pairwise(pw1, confidence=0.95, test.type="dist")

# write.csv(as.data.frame(sum$summary.table), file = here('data/new_TableS5.csv'))




# ------------------------
# trajectory analysis


TA1 <- 
  trajectory.analysis(
    fit1, 
    groups = gdf$taxon, 
    traj.pts = gdf$numerical_stage, 
    print.progress = TRUE)


# shape and magnitude differences between trajectories
summary(TA1, attribute = "SD")
summary(TA1, attribute = "MD")



# ------------------------------
# plot trajectories 



#line colours
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
  plot(TA1, pch = as.numeric(gdf$taxon) +20, 
       bg = traj_col$colour,
       col='white',
       cex=3,
       cex.axis=1.3)
       

#add trajectory lines per species
add.trajectories(
  TP1, 
  traj.pch = c(1,0), 
  traj.lty=c(1,3),
  traj.lwd = c(2,2))

#add legend
legend("topleft",
       levels(gdf$groups), 
       pch =  c(1,0), 
       pt.bg = 1,
       bty ="n",
       cex=1.5)

legend("topright",
       levels(gdf$numerical_stage), 
       pch = 18, 
       title = "Stages",
       col = c("#009E73", "#CC79A7", "#56B4E9", "#E69F00"),
       bty="n",
       pt.cex=3,
       cex=1.5)

#add warp grids
ref <- mshape(proc_data$coords)

par(mar = c(0,0,0,0))
plotRefToTarget(ref, pca_data$pc.shapes$PC1min) #PC1 min
plotRefToTarget(ref, pca_data$pc.shapes$PC1max) #PC1 max
plotRefToTarget(ref, pca_data$pc.shapes$PC2max) #PC2 max
plotRefToTarget(ref, pca_data$pc.shapes$PC2min) #PC2 min






