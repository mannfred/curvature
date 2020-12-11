library(emmeans)
library(here)
library(tidyverse)
here=here::here

# -----------------
# import data ####

# size data
# K_1_2_1 is missing because
# there is no curvature data
size_data <-
  read.csv(
  here("data/epimedium_curv_size_data_nogran.csv"), header=TRUE) %>% 
  tibble() %>% 
  rename(name = species_individual_panicle_flower)
  

# dorsal curvature data
curv_data <- read_rds(path = here('data/RDS_files/curvature_tbl_dorsal.rds'))
  

# check IDs
identical(size_data$name, curv_data$name)



# --------------------
# tidy data


# merge curvature and size data into one tibble
curv_size_data <-
  left_join(curv_data, size_data, by ='name') %>% 
  mutate(species = str_extract(name, "[A-Z]+")) %>%  # extacts uppercase letters from strings
  mutate(new_stage = case_when(
    stage == "E" ~ "C",
    stage == "G" ~ "G",
    stage == "ND" | stage == "F" | stage == "N" | stage == "O" ~ "T",
    stage == "P" | stage == "D" | stage == "A" ~ "A"))

# isolate numeric identifiers
csids <-
  curv_size_data %>% 
  mutate(identity = str_replace(name,  "[a-z]+", "")) %>% #removes alphabet and leaves identifiers 
  mutate(identity = str_replace(identity, substr(identity, 1, 1), ""))  #remove leading "_" from ID strings

# add identifiers to a list
ID_matrix <- str_split_fixed(csids$identity, "_", n=3) #matrix of IDs in SIPC format

# add individual identity
# and convert rad to deg
curv_size_data <-
  csids %>%
  mutate(spp_ind_ID = paste(ID_matrix[,2], species, sep=""))


# -------------------------
# plotting


# colours from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/ 
colour_ids <-
  c(replicate(30, "#E69F00"), #31 E. koreanum samples
    replicate(27, "#56B4E9")) #27 E. violaceum samples

# plot size vs curvature
ggplot(
  data=curv_size_data, 
  aes(x=perimeter, y=total_K)) +
  geom_point(
    aes(colour=factor(species)), size=3) +
  labs(x = "sepal size (mm)", 
       y = "total curvature (degrees)",
       colour = "species") +
  scale_colour_manual(name= "Species", 
                      breaks = c("K", "V"), 
                      labels = c(expression(italic("E. koreanum")), 
                                 expression(italic("E. violaceum"))),
                      values = c("#009E73", "#CC79A7")) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position=c(.85, .85))



# plot stage vs curvature
ggplot() +
  theme(panel.border = element_blank(),    # generally minimal theme
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position=c(.2, .9)) +
  geom_boxplot(
    data=curv_size_data, 
    aes(y=total_K,
        x=factor(new_stage, levels=c("C", "G", "T", "A")),
        fill=species)) +
    labs(x = "developmental stage", 
       y = "total curvature (degrees)") +
    scale_fill_manual(values=c("#009E73", "#CC79A7")) + # fill boxplot colours
    scale_colour_manual(name= "Species", # legend
                    breaks = c("K", "V"), 
                    labels = c(expression(italic("E. koreanum")), 
                               expression(italic("E. violaceum"))))



# ---------------------------
# Fit model for K ~ stage

# random effect (individual) has varies by ~8 std devs
# see ?isSingular for details
model3 <- 
  lmerTest::lmer(total_K ~ new_stage*species + (1|spp_ind_ID), data = curv_size_data)

qqnorm(resid(model3))
qqline(resid(model3))

# pairwise comparisons of means
tukey_results4 <-
  emmeans(model3, pairwise ~ new_stage*species, adjust = "tukey")

# save as data frame
tableS6 <- 
  tukey_results4$contrasts %>% 
  as.data.frame()

 write.csv(tableS6, file=(here("data/new_Table_S6.csv")))


