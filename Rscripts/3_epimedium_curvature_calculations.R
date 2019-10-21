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
        

#inspect raw LMs
str(dorsal_lst)#LMs stored in $coo
dorsal_lst %>% 
  paper %>% 
  draw_curve #draws all curves on one plot (not yet procrustes superimposed)
rapply(dorsal_lst, coo_plot) #recursively plots all curves
coo_plot(dorsal_lst[i]) #for i={1:31} to plot individually

#calculate polynomials
dorsalcurv_lst <- 
  npoly(dorsal_lst$coo, degree=3) #31 polynomial curves estimated 
                                  #"Call:" = function+parameters  used to create the model

#draw polynomials estimates over raw LMs
dorsalcurv_lst[[i]] %>% #for i={1:31} -- LMs must be previously plotted w coo_plot() 
npoly_i() %>%
coo_draw()

#save plot
dev.copy(figure,"Figure_1.png",width=8,height=6,units="in",res=100)
dev.off()


#calculate arc length####

#create polynomial functions from coeffs
coeffs_lst <- 
  sapply(dorsalcurv_lst, function(v) v[1]) 
#extracts the first element ($coeff) from each element ($shp) in the list. 
#sapply stores results as a list of atomic vectors

#function that takes coefficients and constructions functions parameterized by t
aspoly.fun <- 
  function(x) x %>%
  as.numeric() %>%
  as.mpoly() %>% #creates a polynomial object of class 'mpoly'
  print() %>% 
  c("x", .) %>% #parametrize polynomial as t=(x, ax^3+bx^2+cx+d)
  mp() %>% #mpoly
  as.function()

polyfunc_lst <- 
  lapply(coeffs_lst, aspoly.fun) # a list of 31 parameterized polynomial functions


#extracts the lower xy-boundary from b[5], unlists it, and isolates the x-boundary by [1]
#then does the same for b[4] which is the upper xy-boundary stored in dorsal_curv$baseline1
baselines_lst <- 
  dorsalcurv_lst %>% 
  lapply(., 
         function(b) c(unlist(b[5])[1], unlist(b[4])[1])
         )

#calculates arclength for every polynomial bounded by baselines
lengths_lst <- 
  mapply(arclength, polyfunc_lst, 
         baselines_lst %>%
           sapply(., "[[", 1), #min baselines
         baselines_lst %>% 
           sapply(., "[[", 2) #max baselines
         ) %>%
  as.tibble() %>% #row 1 contains arclengths 
  slice(., 1) %>% #keep only row 1
  as.list()





##############################
#calculate total curvature####

poly_lst <- 
  coeffs_lst %>%
  lapply(., polynomial) %>% #convert coeffs to a list of polynomials
  lapply(., as.character) #convert polynomials to character vectors
  

#convert polynomials stored as character strings to quoteless expressions
poly2.fun<- function (p) 
{
  f<- function(x) NULL
  body(f) <- parse(text=p)
  f
}


func_lst<- poly_lst %>%
           lapply(., poly2.fun) #create a list of polynomial functions readable by totalK()


#run trace(fun2form, edit=TRUE) and change width.cutoff to 500L under deparse()

#total curvature function
totalK.fun<-function (x.range, fun) 
{
  stopifnot(is.atomic(x.range))
  if (!is.numeric(x.range)) 
    stop("'x.range' must be a numeric vector!")
  if (length(x.range) != 2) 
    stop("'x.range' must be a vector of length two!")
  if (diff(x.range) < 0) 
    stop("please, reorder 'x.range'.")
  if (!inherits(fun, "function")) 
    stop("'fun' must be a 'function' of x!")
  
  dfun <- deriv3(fun2form(fun), "x", func = TRUE)
  
  if (attr(dfun(x.range[1]), "gradient") == attr(dfun(x.range[2]), #make sure function is not a straight line
                                                 "gradient")) 
    stop("'fun' should not be a linear function of x!") #corrected spelling from "linar"
  x <- seq(x.range[1], x.range[2], length.out = 5000) #splits the x range into 5000 even segments --but does it split the x-axis into 5000 pieces, or the arc length?
  y <- fun(x) 
  gr <- attr(dfun(x), "gradient") #the tangents (first derv) of the 5000 x components, dfun() is defined 7 lines above. The gradient matrix has elements that are the first deriv of a function
  he <- attr(dfun(x), "hessian")[, , "x"] # x is in the third dimension of this object (df?). The hessian matrix has elements that are the second deriv of a function
  k <- abs(he)/(1 + gr^2)^(3/2) #there are possibly 5000 curvature points stored in k... could sum them? also, see: https://en.wikipedia.org/wiki/Curvature "curvature of the graph of a function"
}


#calculate curvature many times along many curves
curvature_tbl <- 
  totalK.fun %>%
  mapply(., baselines_lst, func_lst) %>%
  as.tibble() %>%
  rapply(., sum) %>% #integration of f dx 
  as.tibble() %>%
  gather()


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
  mutate(adjusted_curvature = total_curvature/arclength) %>% #new column is adjusted curvature
 # write.csv(., here("data/epimedium_curvature_koreanum.csv")) 
  





#procrustes superimposition####

dorsal_procrustes<-fgProcrustes(dorsal) 
dorsal_procrustes%>% paper %>% draw_curve #draws all procrustes-superimposed curves
#two are inverted?
rapply(dorsal_procrustes, coo_plot) 


#fit polynomials
curves<-opoly(dorsal_procrustes) #degree=5 by default
