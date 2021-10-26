# import PCs of shape space (see: "3_geomorph_analysis.R")
# to regress against total curvature estimates

library(tidyverse)
library(effectsize)
library(here)
library(geomorph)
library(lme4)
library(lmerTest)
library(vegan)
here=here::here

# ----------------------------------------
# import procrustes alignment and PCA data

proc_data <- readRDS(file=(here("data/derived_data/RDS_files/proc_data.rds")))

pca_data <- readRDS(file=here("data/derived_data/RDS_files/pca_data.rds"))


# ---------------------------------------
# reproducing plotTangentSpace() algorithm
# to recreate PCA shape space
# p = number of landmarks
# k = number of dimensions (2)
# n = number of specimens

# convert 3D npk array into 2D n[pk] array 
# shape_vars is the shape data with which pca is implemented
# this was checked by modifying plotTangentSpace() to take
# shape_vars: the PCA results are identical
shape_vars <- two.d.array(proc_data$coords)



# tol value calculated by modifying plotTangentSpace()
# to include print(tol) which gives 0.005 for this dataset
pc_res <- prcomp(shape_vars, tol = 0.005) 



# check to see that shape_vars is equivalent
# dataset that plotTangentSpace() processes
identical(pc_res$x, pca_data$pc.scores) #TRUE

pc_res_summary <- 
  summary(pc_res)$importance %>% 
  as.data.frame()


# colours from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/ 
colour_ids <-
  c(replicate(30, "#E69F00"), #31 E. koreanum samples
    replicate(27, "#56B4E9")) #27 E. violaceum samples

# plot PCA shape space
plot(
  pc_res$x, 
  type= 'points', 
  col = colour_ids, 
  pch = 21, 
  bg = colour_ids,
  xlab = 
    paste0('PC1', ' ', '(', round(100*pc_res_summary$PC1[2], digits=2), '%)' ),
  ylab = 
    paste0('PC2', ' ', '(', round(100*pc_res_summary$PC2[2], digits=2), '%)' ),
  text=title(main="Unconstrained PCA of shape variables"))


# compare to plotTangentSpace()
plot(pca_data$pc.scores[,1:2], 
     col = colour_ids,
     pch=21,
     bg=colour_ids) 



# -------------------
#redundancy analysis

#if two.d.array() outputs are 'shape variables',
#then rda(shape_vars ~ adj_curv) is a
#redundancy analysis of shape constrained by curvature

#import dorsal curvature data

curv_data <- 
  read_rds(path = here('data/derived_data/RDS_files/spline_curvature_tbl_dorsal.rds'))

 

#TRUE = all sample IDs match between shape matrix and
#curvature vector
identical(as.character(curv_data$name), unlist(dimnames(shape_vars)[1])) 

#redundancy analysis
const_ord <- 
  rda(shape_vars ~ curv_data$total_K)

const_ord_summary <- 
  summary(const_ord)$cont$importance %>% 
  as.data.frame()

#plotting
plot(
  const_ord, 
  type = 'points', 
  col = colour_ids, 
  pch =21, 
  bg =colour_ids,
  xlab = paste0(
    'RDA1: total curvature', ' ', '(', 
    round(100*const_ord_summary$RDA1[2], digits=2), 
    '%)' ),
  ylab = paste0(
    'PC1', ' ', '(', 
    round(100*const_ord_summary$PC1[2], digits=2), 
    '%)' ))

title(main= "RDA with single explantory var.")

vegan::anova.cca(const_ord, by='axis')



# ------------------------------

# inspecting and comparing individual
# ordination axes


rda_scores <- scores(const_ord, choices=c(1:4))$sites 

#plot PC1(rda) vs PC2(rda)
plot(
  rda_scores[,2], 
  rda_scores[,3],
  type = 'points', 
  xlab = paste0(
    'PC1(r)', ' ', '(', 
    round(100*const_ord_summary$PC1[2], digits=2), 
    '%)' ),
  ylab = paste0(
    'PC2(r)', ' ', '(', 
    round(100*const_ord_summary$PC2[2], digits=2), 
    '%)' ),
  col = colour_ids, 
  pch=21, 
  bg=colour_ids)

title(main= "PCA of Y ordinated by curvature")


# ---------------------
# partial PCA 


# pPCA equivalent to rda(shape_vars ~ shape_vars + adj_crv)
shape_residuals <- residuals(lm(shape_vars ~ curv_data$total_curvature))

pPCA <- prcomp(shape_residuals)

# colours from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/ 
colour_ids <-
  c(replicate(30, "#E69F00"), #31 E. koreanum samples
    replicate(27, "#56B4E9")) #27 E. violaceum samples

plot(pPCA$x, col = colour_ids, pch=21, bg=colour_ids)

pPCA_summary <- summary(pPCA)

arrayspecs(pPCA$, 31, 2) 

plotRefToTarget(ref, pPCA$x) #PC1 min

# plot RDA1 vs pPC1 
plot(
  rda_scores[,1], 
  pPCA$x[,1], 
  col = colour_ids, 
  pch=21, 
  bg=colour_ids)



# --------------------------
#partial PCA of shape ~ RDA1


shape_res2 <- residuals(lm(shape_vars ~ rda_scores[,1]))

pPCA2 <- prcomp(shape_res2)

pPCA2_summary <- 
  summary(pPCA2)$importance %>% 
  as.data.frame()

plot(
  rda_scores[,1], 
  pPCA2$x[,1], 
  col = colour_ids, 
  pch=21, 
  bg=colour_ids,
  xlab= paste0(
    'RDA1', ' ', '(', 
    round(100*const_ord_summary$RDA1[2], digits=2), 
    '%)' ),
  ylab = paste0(
    'pPC1', ' ', '(', 
    round(100*pPCA2_summary$PC1[2], digits=2), 
    '%)' ))

title(main= "RDA1 vs pPC1")



# -------------------------------
# are PC1 of shape space and 
# adjusted curvature correlated?


# isolate numeric identifiers
csids <-
  curv_data %>% 
  mutate(identity = str_replace(name,  "[a-z]+", "")) %>% #removes alphabet and leaves identifiers 
  mutate(identity = str_replace(identity, substr(identity, 1, 1), ""))  #remove leading "_" from ID strings

# add identifiers to a list
ID_matrix <- str_split_fixed(csids$identity, "_", n=3) #matrix of IDs in SIPC format

# create group vector
group_ids <-
  c(replicate(30, "koreanum"),
    replicate(27, "violaceum")) %>% 
  factor()

# put species epithets back in
curv_data <-
  curv_data %>% 
  mutate(taxon = group_ids)

# add individual identity
curv_data <-
  csids %>%
  mutate(spp_ind_ID = paste(curv_data$taxon, ID_matrix[,2], sep=""))



# set model variables
PC1 <- pca_data$pc.scores[,1]
PC2 <- pca_data$pc.scores[,2]
PC3 <- pca_data$pc.scores[,3]
dors_curv <- curv_data$total_K
indiv <- curv_data$spp_ind_ID

# models with random effect for indiv
model1 <- 
  lmerTest::lmer(PC1 ~ dors_curv*group_ids + (1|indiv))

anova(model1) #F=8.8303, NumDF = 1, DenDF = 52.604
effectsize::F_to_eta2(20.1468, 1, 52.657) #eta2 = 0.14



model2 <-
  lmerTest::lmer(PC2 ~ dors_curv*group_ids + (1|indiv))

anova(model2)
effectsize::F_to_eta2(60.814, 1, 52.999) #eta2=0.69




# ----------------------------
# plot total curvature vs PC2

mydata <- tibble(dors_curv, PC2)

ggplot(data = mydata, aes(x = dors_curv, y = PC2, colour=colour_ids)) +
  geom_point(size=4) +
  geom_smooth(method = 'lm', se = FALSE, size=1.5) +
  labs(
    y = "PC2 of Shape Variation (25.6%)",
    x = "total curvature (degrees)") +
  scale_colour_manual(
    name = "Species", 
    breaks = c("#E69F00", "#56B4E9"),
    labels = c(expression(italic("E. koreanum")), expression(italic("E. violaceum"))),
    values = c("#009E73", "#CC79A7")) +
  theme_classic() + 
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=14)) 

