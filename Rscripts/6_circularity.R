library(here)
library(tidyverse)



# -----------------
# import data
circularity <- 
  readRDS(file=here('data/RDS_files/circularity.rds')) %>%  
  tibble() %>% 
  rename(circularity = ".")

size_data<-
  read.csv(
    here("data/epimedium_curv_size_data_nogran.csv"), header=TRUE) %>% 
  tibble() %>% 
  rename(name = species_individual_panicle_flower) %>% 
  mutate(new_stage = case_when(
    stage == "E" ~ "C",
    stage == "G" ~ "G",
    stage == "ND" | stage == "F" | stage == "N" | stage == "O" ~ "T",
    stage == "P" | stage == "D" | stage == "A" ~ "A"))

# dorsal curvature data
curv_data <-
  read.csv(
    here("data/epimedium_adj_curvature_dorsal.csv"), header=TRUE) %>% 
  tibble() %>% 
  select(2:4) 

curv_circ_data <- cbind(circularity, curv_data, size_data[,-1])

# ------------------------------
# plot stage vs circularity

ggplot(
  data=curv_circ_data, 
  aes(x=new_stage, y=circularity)) +
  geom_boxplot(
    aes(fill = factor(species), #colour by species
        factor(new_stage, levels=c("C", "G", "T", "A")))) + #reorder
  labs(x = "developmental stage", 
       y = "circularity (unit-less)",
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
        legend.position=c(.2, .2))


# -----------------------
# fit model

library(emmeans)

# isolate numeric identifiers
csids <-
  curv_circ_data %>% 
  mutate(identity = str_replace(name,  "[a-z]+", "")) %>% #removes alphabet and leaves identifiers 
  mutate(identity = str_replace(identity, substr(identity, 1, 1), ""))  #remove leading "_" from ID strings

# add identifiers to a list
ID_matrix <- str_split_fixed(csids$identity, "_", n=3) #matrix of IDs in SIPC format

# add individual identity
curv_circ_data <-
  csids %>%
  mutate(spp_ind_ID = paste(ID_matrix[,2], species, sep=""))


model4 <- 
  lmerTest::lmer(circularity ~ new_stage*species + (1|spp_ind_ID), data = curv_circ_data)

qqnorm(resid(model4))
qqline(resid(model4))


tukey_results4 <-
  emmeans(model4, pairwise ~ new_stage*species, adjust = "tukey")

write.csv(tukey_results4$contrasts, file=here('data/new_TableS8.csv'))




# ---------------------------------
# look for local curvature anomalies

lowk <- total_curvature(baselines_list[[13]], poly_list[[13]], 500)
hihk <- total_curvature(baselines_list[[7]], poly_list[[7]], 500)

lowf <- curvr::as_function(poly_list[[13]])
hihf <- curvr::as_function(poly_list[[7]])

# ------

x <- seq(from = -2.6439, to = 1.4663, by = 0.1)
y <- lowf(x)

plot.new()
plot.window(xlim=c(-2.6439, 1.4663), ylim =c(-1.5, 1.5))
points(x, y)

plot(1:501, lowk$k)
plot(x,y)

# ---------
x2 <- seq(from = -4.376075, to = 2.595412, by = 0.1)
y2 <- hihf(x2)

plot(x2, y2)
plot(1:501, hihk$k)


plot.new()
plot.window(xlim=c(-4.37, 2.595), ylim =c(-4, 4))
points(x2, y2)

