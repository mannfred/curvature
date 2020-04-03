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
  select(2:5)

# merge curvature data into one tibble
curv_size_data <-
  left_join(curv_data, size_data, by ='name') %>% 
  as_tibble() %>% 
  mutate(species = str_extract(name, "[A-Z]+")) # extacts uppercase letters from strings
  

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
  mutate(indiv_ID = ID_matrix[,2]) #add individual identity



# plot sepal size vs curvature
ggplot(
  data=curv_size_data, 
  aes(x=sepal_size_mm, y=log(adjusted_curvature))) +
  geom_point(
    aes(colour=factor(species)), size=3) +
  labs(x = "Sepal length (mm)", 
       y = "log(total mean curvature) (degrees/mm)",
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
  


##################################
#Fit model for log(K) ~ sepal size


model3 <- 
  lmerTest::lmer(log(adjusted_curvature) ~ Species_epithet + new_stage + (1|indiv_ID), data = data5)

qqnorm(resid(model1))
qqline(resid(model1))

tukey_results<-
  emmeans(model1, list(pairwise ~ Species_epithet + new_stage), adjust = "tukey")



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
