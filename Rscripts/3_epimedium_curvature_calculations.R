########################################
#measure curvature from epimedium photos


library(here)
library(Momit)
library(Momocs)
library(geomorph)
library(tidyverse)
library(mpoly)
library(pracma)
library(soilphysics)
library(polynom)

#import dorsal data####
dorsal_lst <- 
  here("data/epimedium_photos/koreanum/koreanum_dorsal_appended.TPS") %>% #mounts tps file
  readland.tps(specID='imageID') %>% #geomorph func, auto-applies the scaling factor to LMs
  a2l() %>% # Momit func, converts array to list
  Ldk() #adds "$shp{i}" as list element headers
        


#alternative to above
#import data
# dorsal_data<-
#   here("data/epimedium_photos/koreanum/koreanum_dorsal_appended.TPS") %>%
#   import_tps() 

#not functional yet..
# #format data and apply scaling
# dorsal_lst<-
#   map2(dorsal_data$coo, dorsal_data$scale, function(x, y) x*y) %>%   #apply scaling
  

#inspect raw LMs
# str(dorsal_lst)#LMs stored in $coo
# dorsal_lst[1] %>%
#   paper %>%
#   draw_curve #draws all curves on one plot (not yet procrustes superimposed)
# rapply(dorsal_lst, coo_plot) #recursively plots all curves
# coo_plot(dorsal_lst[1]) #for i={1:31} to plot individually

#calculate polynomials from Momocs::npoly
dorsalcurv_lst <- 
  npoly(dorsal_lst$coo, degree=2) #Momocs::npoly object with coeffs, baselines, etc
                                  #"Call:" = function+parameters  used to create the model








#draw polynomials estimates over raw LMs
# dorsalcurv_lst[[1]] %>% #for i={1:31} -- LMs must be previously plotted w coo_plot() 
# npoly_i() %>%
# coo_draw()
# 
# #save plot
# dev.copy(figure,"Figure_1.png",width=8,height=6,units="in",res=100)
# dev.off()




#########################
#calculate arc length####

#use param()
paramfun_lst <- 
  lapply(dorsalcurv_lst, param) # a list of 31 parameterized polynomial functions


#extracts the lower xy-boundary from b[5], unlists it, and isolates the x-boundary by [1]
#then does the same for b[4] which is the upper xy-boundary stored in dorsal_curv$baseline1
baselines_lst <- 
  dorsalcurv_lst %>% 
  lapply(., 
         function(b) c(unlist(b[5])[1], unlist(b[4])[1])
         )

#calculates arclength for every polynomial bounded by baselines
lengths_lst <- 
  mapply(arclength, paramfun_lst, #applies arclength() to all elements in polyfunc_lst
         baselines_lst %>%
           sapply(., "[[", 1), #min baselines, "[[" is a subsetting function
         baselines_lst %>% 
           sapply(., "[[", 2) #max baselines
         ) %>%
  as_tibble() %>% #row 1 contains arclengths 
  slice(., 1) %>% #keep only row 1
  as.list()





##############################
#calculate total curvature####

#add assign("dx_list", dX, envir = .GlobalEnv) to trace(arclength, edit=TRUE) to inspect objects created during arclength() calculations

#calculate curvature many times along many curves
curvature_tbl <- 
  totalK_fun %>%
  mapply(., baselines_lst[[1]], func_lst[[1]], dorsalcurv_lst[[1]]) %>%
  as.tibble() %>%
  rapply(., sum) %>% #integration of f dx 
  as.tibble() %>%
  gather()

saveRDS(curvature_tbl, file=here("data/RDS_files/curvature_tbl.rds"))

#merge curvature and arclength information
alltogether_tbl <- 
  lengths_lst %>% #list of arclengths
  unlist() %>%
  as.tibble() %>%
  gather() %>% #unpivot column names to row names
  {bind_cols(                 #bind arclength column and curvature column
    dplyr::select(., value ), #mask Momocs::select
    dplyr::select(curvature_tbl, value),
    dplyr::select(  
              read.csv(
                  here("data/epimedium_curv_size_data.csv"), 
                  header=TRUE) %>% #fetch ID tags from data.csv
              dplyr::select(species_individual_panicle_flower) %>% #isolate ID column
              slice(., 20:50) %>%   #20:50 for koreanum, 51:77 for violaceum, 1:19 for grandiflorum
              as.tibble(),
                  species_individual_panicle_flower
                  )
            )
  } %>% 
  rename(arclength=value, 
         total_curvature=value1) %>% #old colnames were "value" and "value1"
  mutate(adjusted_curvature = total_curvature/arclength) #%>% #new column is adjusted curvature
 # write.csv(., here("data/epimedium_curvature_koreanum.csv")) 
  





#procrustes superimposition####

dorsal_procrustes<-fgProcrustes(dorsal) 
dorsal_procrustes%>% paper %>% draw_curve #draws all procrustes-superimposed curves
#two are inverted?
rapply(dorsal_procrustes, coo_plot) 


#fit polynomials
curves<-opoly(dorsal_procrustes) #degree=5 by default
