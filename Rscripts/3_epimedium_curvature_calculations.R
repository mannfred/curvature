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
  

# coeffs_lst <- 
#   sapply(dorsalcurv_lst, function(v) v[1]) 
    
#written by user: 李哲源 at https://stackoverflow.com/questions/40438195/function-for-polynomials-of-arbitrary-order-symbolic-method-preferred 
express <- function (npoly, expr = TRUE) {
  coeffs <- npoly[[1]] #extract coefficients from npoly list
  stringexpr <- paste("x", seq_along(coeffs) - 1, sep = " ^ ")
  stringexpr <- paste(stringexpr, coeffs, sep = " * ")
  stringexpr <- paste(stringexpr, collapse = " + ")
  if (expr) return(parse(text = stringexpr))
  else return(stringexpr)
}

#exp_lst<-lapply(dorsalcurv_lst, express) corresponding lapply() function for express()
#how to create a fuction programatically
#https://stackoverflow.com/questions/12982528/how-to-create-an-r-function-programmatically
#https://stackoverflow.com/questions/9345373/as-alist-character
#should look at https://github.com/dkahle/mpoly/tree/master/R 


param <- function(npoly, expr = TRUE) {
  coeffs <- npoly[[1]] #extract coefficients from npoly list
  stringchar <- paste("x", seq_along(coeffs) - 1, sep = " ^ ")
  stringchar <- paste(stringchar, coeffs, sep = " * ")
  stringchar <- paste(stringchar, collapse = " + ")
  paramchar <- paste("x", stringchar, sep = " , " )
  bodyexp <- eval(parse(text = paste("alist(", paramchar, ")")))
  
  f<- function(x) NULL
  body(f) <- paste(bodyexp[1], bodyexp[2], sep = ",")
  return(f) 
}
  

  




#  deriv_lst<-lapply(exp_lst, deriv3, "x")




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

paramfun_lst <- 
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
  mapply(arclength, polyfunc_lst, #applies arclength() to all elements in polyfunc_lst
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

poly_lst <- 
  coeffs_lst %>%
  lapply(., polynomial) %>% #convert coeffs to a list of polynomials (character strings, not functions)
  lapply(., as.character) #convert polynomials to character vectors
  

#convert polynomials stored as character strings to quoteless expressions
poly2_fun<- function (p) 
{
  f<- function(x) NULL
  body(f) <- parse(text=p)
  f
}


func_lst<- poly_lst %>%
           lapply(., poly2_fun) #create a list of polynomial functions readable by totalK()


#run trace(fun2form, edit=TRUE) and change width.cutoff to 500L under deparse()

#see: https://github.com/mannfred/curvature/blob/8226cd29b88fa83b2c9ce52f16058c73c4921865/Rscripts/3_epimedium_curvature_calculations.R
#for un-screwy version
#total curvature function
totalK_fun<-function (x_range, poly_list, subdiv) 
{
  stopifnot(is.atomic(x_range)) # is.atomic checks that x.range cannot be a list or expression
  if (!is.numeric(x_range)) 
    stop("'x_range' must be a numeric vector!")
  if (length(x_range) != 2) 
    stop("'x_range' must be a vector of length two!")
  if (diff(x_range) < 0) #calculates the difference between x1 and x2 to ensure it's >0
    stop("please reorder 'x_range'.")
  # if (!inherits(fun, "function")) 
  #   stop("'fun' must be a 'function' of x!")
  
  # dfun <- deriv3(fun2form(fun), "x", func = TRUE) #func=TRUE returns a function
  
  exp_list <- lapply(poly_list, f) #a list of polynomial expressions to be read by deriv3()
  dfun <- deriv3(exp_list, "x", func=TRUE) %>% unname()
  
  if (attr(dfun(x_range[1]), "gradient") == attr(dfun(x_range[2]), "gradient"))  #make sure function is not a straight line
    stop("'fun' should not be a linear function of x!") 
  
  iter<- seq(0, 1, by=1/subdiv) #create vector of subdivisions to calculate arclength parameter
  arcfun_lst<- list() #empty bin
  b<- arclength(param_fun, x_range[1], x_range[2])$length #arc length of t-parameterized function
  
  for(i in seq_along(iter)){ 
    arcfun_lst[[i]] <- 
      local({
        b_sub<-iter[i]*b
        function(u) arclength(param_fun, x_range[1], u)$length - b_sub
      }) 
  }
  
  root_find<- function(x) uniroot(x, x_range)$root #root-finding function
  
  x <- sapply(arcfun_lst, root_find) #find roots for a list of 
  
  poly2_fun<- function (p) #function for converting expressions to functions
  {
    f<- function(x) NULL
    body(f) <- parse(text=p)
    f
  } 
  
  func_lst<- exp_list %>%
    lapply(., poly2_fun) #create a list of polynomial functions readable by totalK()
  
  y <- func_lst(x) 
  gr <- attr(dfun(x), "gradient") #the tangents (first derv) of the x_n components, dfun() is defined above. The gradient matrix has elements that are the first deriv of a function
  he <- attr(dfun(x), "hessian")[, , "x"] # x is in the third dimension of this object (df?). The hessian matrix has elements that are the second deriv of a function
  k <- abs(he)/(1 + gr^2)^(3/2) #has n=subdiv measurements of k
  k_total<-(sum(k) /subdiv) *(180/pi) #add all measurements of k, rescale depending on #of subdivisions, and convert from rad to degrees
}




#add assign("dx_list", dX, envir = .GlobalEnv) to trace(arclength, edit=TRUE) to inspect objects created during arclength() calculations
# b<- arclength(param_fun, x_range[1], x_range[2])$length


############testing START######

#in totalK_fun: can be simplified if I can find out how to convert functions into parameterized functions

x<-seq(0, 0.9999, by=0.0001) 


#circle 

circlefun<-function(x) (1-(x^2))^0.5

dfun2<- deriv3(fun2form(function(x) (1-(x^2))^0.5), "x", func=TRUE)
param_fun<-function(t) c((1-(t^2))^0.5, t)
fun<-function(x) (1-(x^2))^0.5

gr2 <- attr(dfun2(x), "gradient") #computes 1st derivative bw x=0 to x=1

he2 <- attr(dfun2(x), "hessian")[ , , "x"] #computes 2nd derivative bw x=0 to x=1


k2 <- abs(he2)/(1 + gr2^2)^(3/2)
(sum(k2) * 0.0001) #*(180/pi)




############testing END######


#calculate curvature many times along many curves
curvature_tbl <- 
  totalK_fun %>%
  mapply(., baselines_lst[[1]], func_lst[[1]], polyfunc_lst[[1]]) %>%
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
