library(here)
library(tidyverse)
here=here::here
#for E. grandiflorum only (all 3 spp. done simultaneously below)####

#combine alignment matrix with size data 

#alignment data
stages_days_subset <- #only import E. grandiflorum data
  readRDS("data\stages_days_sort_grandiflorum.rds") %>%
  na_if("-") %>% 
  drop_na() %>% #remove rows with "-" 
  group_by(ID) %>% #create groups by flower
  mutate(row_ID = row_number()) %>%
  unite(unique_ID,  ID, row_ID, sep="_") %>% #create unique_ID column
  ungroup()

#size data
data5 <- readRDS("epimedium_growth_data_pivot_redefined_stages.rds")

data5_subset<- data5 %>%
  filter(Species_epithet == 'grandiflorum') %>%
  select(c(1, size)) %>%
  group_by(Species_Individual_Panicle_Flower) %>%
  mutate(row_ID = row_number()) %>%
  unite(unique_ID,  Species_Individual_Panicle_Flower, row_ID, sep="_") %>% #create unique_ID column
  ungroup()

joiner_grandiflorum<-
  stages_days_subset %>%
  left_join( data5_subset, by = "unique_ID")

#plot single species
ggplot(
  data=joiner_grandiflorum, 
  aes(x=elapsed_days, y=size)
) +
  geom_point(
    aes(x=elapsed_days, y=size))




##########################################
# make scatterplot with all 3 species ####

#import elapsed days data
gran_stages_days <- 
  readRDS("stages_days_sort_grandiflorum.rds") %>%
  na_if("-") %>%
  drop_na() %>% #remove rows with "-"
  group_by(ID) %>% #create groups by flower
  mutate(row_ID = row_number()) %>%
  unite(unique_ID,  ID, row_ID, sep="_") %>% #create unique_ID column
  ungroup()

kore_stages_days <-
  readRDS(here("data/RDS_files/stages_days_sort_koreanum.rds")) %>%
  na_if("-") %>%
  drop_na() %>% #remove rows with "-"
  group_by(ID) %>% #create groups by flower
  mutate(row_ID = row_number()) %>%
  unite(unique_ID,  ID, row_ID, sep="_") %>% #create unique_ID column
  ungroup()

viol_stages_days <-
  readRDS(here("data/RDS_files/stages_days_sort_violaceum.rds")) %>%
  na_if("-") %>%
  drop_na() %>% #remove rows with "-"
  group_by(ID) %>% #create groups by flower
  mutate(row_ID = row_number()) %>%
  unite(unique_ID,  ID, row_ID, sep="_") %>% #create unique_ID column
  ungroup()

#import size data
data5 <- readRDS(here("data/RDS_files/epimedium_growth_data_pivot_redefined_stages.rds"))

size_data <- 
  data5 %>%
  filter(str_detect(Species_Individual_Panicle_Flower, "[kv]")) %>% #excludes grandiflorum
  select(c(1, size)) %>% #choose columns with species IDs and size data
  group_by(Species_Individual_Panicle_Flower) %>%
  mutate(row_ID = row_number()) %>% #create unique identifiers by row number
  unite(unique_ID,  Species_Individual_Panicle_Flower, row_ID, sep="_") %>% #create unique_ID column
  ungroup()


#join E. koreanum and E. violaceum
joiner <- 
  list(size_data, kore_stages_days, viol_stages_days) %>%
  reduce(full_join, by = "unique_ID") %>%
  mutate(days_elapsed = coalesce(elapsed_days.x, elapsed_days.y)) %>% #create merged days column
  mutate(flower_stage = coalesce(stage.x, stage.y)) %>% #create merged stage column
  select(-c(3:6)) %>% #remove old columns
  mutate(species_epithet = str_extract(unique_ID, "[a-z]+")) %>%  #extacts lowercase letters from strings
  group_by(species_epithet)

#plot
ggplot(
  data=joiner, 
  aes(x=days_elapsed, y=size)) +
geom_point(
  aes(colour=factor(species_epithet)), size=3) +
labs(x = "Elapsed time (days)", 
     y = "Distance between sepals (mm)") +
scale_colour_manual(name= "Species", 
                    breaks = c("koreanum", "violaceum"), 
                    labels = c(expression(italic("E. koreanum")), 
                               expression(italic("E. violaceum"))),
                    values = c("#009E73", "#CC79A7")) +
theme_bw() + 
theme(panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"), 
      legend.position=c(.25, .65),
      legend.text = element_text(size=14),
      legend.title=element_text(size=14))





# -----------------------------
# fit models for stage ~ days

#isolate numeric identifiers
sppids <-
  joiner %>% 
  mutate(identity = str_replace(unique_ID,  "[a-z]+", "")) %>% #removes alphabet and leaves identifiers 
  mutate(identity = str_replace(identity, substr(identity, 1, 1), ""))  #remove leading "_" from ID strings

#add identifiers to a list
ID_matrix <- str_split_fixed(sppids$identity, "_", n=3) #matrix of IDs in SIPC format

#add identifiers back to tibble 
joiner <-
  sppids %>%
  ungroup() %>% 
  mutate(spp_ind_ID = paste(ID_matrix[,1], species_epithet, sep='')) #add indiv_ID

# mixed effects model
model2 <- 
  lmerTest::lmer(days_elapsed ~ flower_stage*species_epithet + (1|spp_ind_ID), data = joiner)

qqnorm(resid(model2))
qqline(resid(model2))

tukey_results2<-
  emmeans(model2, list(pairwise ~ species_epithet*flower_stage), adjust = "tukey")

# ----------------------
# plot stage vs days

emmip(model2, species_epithet~flower_stage, type="response") +
  geom_point(
    aes(x = flower_stage, y = days_elapsed, colour = species_epithet), 
    data = joiner, 
    pch = 20, 
    alpha=0.5, 
    size=5, 
    position=position_jitterdodge(dodge.width=0.3)) +
  labs(y = "Days elapsed",
       x = "Developmental stage") +
  scale_x_discrete(limits=c("C", "G", "T", "A")) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.20)) +
  stat_summary(
    fun.y = mean, 
    geom = "errorbar", 
    aes(ymax = ..y.., ymin = ..y.., group = factor(species_epithet)),
    width = 0.5, 
    linetype = "solid", 
    position = position_dodge(),
    size=1.5,
    alpha=0.65) +
  scale_colour_manual(name= "Species", 
                      breaks = c("koreanum", "violaceum"), 
                      labels = c(expression(italic("E. koreanum")), 
                                 expression(italic("E. violaceum"))),
                      values = c("#009E73", "#CC79A7")) 



# ---------------------------
# model stage vs size


# mixed effects model
model1 <- 
  lmerTest::lmer(size ~ flower_stage*species_epithet + (1|spp_ind_ID), data = joiner)

qqnorm(resid(model1))
qqline(resid(model1))

tukey_results1<-
  emmeans(model1, list(pairwise ~ species_epithet*flower_stage), adjust = "tukey")


#plot model

emmip(model1, species_epithet~flower_stage, type="response") +
  geom_point(
    aes(x = flower_stage, y = size, colour = species_epithet), 
    data = joiner, 
    pch = 20, 
    alpha=0.5, 
    size=5, 
    position=position_jitterdodge(dodge.width=0.3)) +
  labs(y = "Sepal size (mm)",
       x = "Developmental stage") +
  scale_x_discrete(limits=c("C", "G", "T", "A")) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.20)) +
  stat_summary(
    fun.y = mean, 
    geom = "errorbar", 
    aes(ymax = ..y.., ymin = ..y.., group = factor(species_epithet)),
    width = 0.5, 
    linetype = "solid", 
    position = position_dodge(),
    size=1.5,
    alpha=0.65) +
  scale_colour_manual(name= "Species", 
                      breaks = c("koreanum", "violaceum"), 
                      labels = c(expression(italic("E. koreanum")), 
                                 expression(italic("E. violaceum"))),
                      values = c("#009E73", "#CC79A7")) 
