library(emmeans)
library(here)
library(tidyverse)



#import data ####

# size data
size_data<-
  read.csv(
  here("data/epimedium_curv_size_data.csv"), header=TRUE) %>% 
  slice(20:77) %>% 
  rename(name = species_individual_panicle_flower)
  

# curvature data
curv_data <-
  read.csv(
  here("data/epimedium_adj_curvature.csv"), header=TRUE) %>% 
  select(2:5) #drop column "X"
 

# merge curvature data into one tibble
curv_size_data <-
  left_join(curv_data, size_data, by ='name') %>% 
  as_tibble() %>% 
  mutate(species = str_extract(name, "[A-Z]+")) %>%  # extacts uppercase letters from strings
  mutate(new_stage = case_when(
    stage == "E" ~ "C",
    stage == "G" ~ "G",
    stage == "ND" | stage == "F" | stage == "N" | stage == "O" ~ "T",
    stage == "P" | stage == "D" | stage == "A" ~ "A"))

#isolate numeric identifiers
csids <-
  curv_size_data %>% 
  mutate(identity = str_replace(name,  "[a-z]+", "")) %>% #removes alphabet and leaves identifiers 
  mutate(identity = str_replace(identity, substr(identity, 1, 1), ""))  #remove leading "_" from ID strings

#add identifiers to a list
ID_matrix <- str_split_fixed(csids$identity, "_", n=3) #matrix of IDs in SIPC format

#add identifiers back to tibble 
curv_size_data <-
  csids %>%
  mutate(spp_ind_ID = paste(ID_matrix[,2], species, sep="")) #add individual identity



# plot sepal size vs curvature
ggplot(
  data=curv_size_data, 
  aes(x=sepal_size_mm, y=adjusted_curvature)) +
  geom_point(
    aes(colour=factor(species)), size=3) +
  labs(x = "Sepal length (mm)", 
       y = "total mean curvature (degrees/mm)",
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
        legend.position=c(.85, .65))
  


#############################
#Fit model for K ~ sepal size


model3 <- 
  lmerTest::lmer(adjusted_curvature ~ new_stage*species + (1|spp_ind_ID), data = curv_size_data)

qqnorm(resid(model3))
qqline(resid(model3))


emmeans(model3, list(pairwise ~ new_stage*species), adjust = "tukey")





##################################
#Fit model for K ~ centroid size

#import geomorph data frame w centroid info
gdf <- readRDS(file=here("data/RDS_files/geomorph_data_frame.rds"))

#add to main tibble
curv_size_data <-
  curv_size_data %>% 
  mutate(Csize = gdf$Csize)

#
model4 <- 
  lmerTest::lmer(adjusted_curvature ~ Csize*species + (1|spp_ind_ID), data = curv_size_data)

qqnorm(resid(model4))
qqline(resid(model4))


emmeans(model4, list(pairwise ~ Csize*species), adjust = "tukey")

#very similar trend to sepal size vs K
#therefore, sepal size is an appropriate proxy for centroid size



##############################
# plot sepal size vs arclength
ggplot(
  data=curv_size_data, 
  aes(x=sepal_size_mm, y=arclength)) +
geom_point(
  aes(colour=factor(species)), size=3) +
labs(x = "Sepal length (mm)", 
     y = "Total mean curvature (degrees/mm)",
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
      legend.position=c(.85, .65))
               





# comapre arc lengths at maturity (T stage and beyond)
mature_data <-
  curv_size_data %>% 
  filter(sepal_size_mm >= 23.0) %>% #22.9 mm is the upper CI for "T" stage
  group_by(species)

ggplot(mature_data, aes(species, arclength)) +
  geom_boxplot()

t.test(arclength ~ species, data = mature_data)
